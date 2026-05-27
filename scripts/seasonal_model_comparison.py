from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib-cache"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from netCDF4 import Dataset, num2date
import numpy as np
import pandas as pd
import xarray as xr


PROCESSED = ROOT / "data" / "processed"


@dataclass(frozen=True)
class Season:
    name: str
    label: str
    start: str
    end: str


DEFAULT_SEASONS = [
    Season("early", "Early season", "2024-09-01", "2024-12-31"),
    Season("midwinter", "Mid-winter", "2025-01-01", "2025-03-31"),
    Season("melt", "Melt season", "2025-04-01", "2025-06-01"),
    Season("full", "Full season", "2024-09-01", "2025-06-01"),
]


def load_snow_depths(processed_dir: Path = PROCESSED) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
    """Load snow depth from all three products and convert to metres."""
    fsm_ds = xr.open_dataset(processed_dir / "fsm2_sd_all.nc")
    senorge_ds = xr.open_dataset(processed_dir / "senorge_all.nc")
    snowpack_ds = xr.open_dataset(processed_dir / "snowpack_all.nc")

    fsm_sd = fsm_ds["snow_depth"] / 1000.0
    senorge_sd = senorge_ds["snow_depth"] / 100.0
    snowpack_sd = snowpack_ds["HS_mod"] / 100.0

    fsm_sd.attrs["units"] = "m"
    senorge_sd.attrs["units"] = "m"
    snowpack_sd.attrs["units"] = "m"

    return snowpack_sd, fsm_sd, senorge_sd


def load_snowpack_depth(processed_dir: Path = PROCESSED) -> xr.DataArray:
    """Load SNOWPACK snow depth and convert to metres."""
    snowpack_ds = xr.open_dataset(processed_dir / "snowpack_all.nc")
    snowpack_sd = snowpack_ds["HS_mod"] / 100.0
    snowpack_sd.attrs["units"] = "m"
    return snowpack_sd


def prepare_station_series(
    snowpack_sd: xr.DataArray,
    fsm_sd: xr.DataArray,
    senorge_sd: xr.DataArray,
) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
    """Resample SNOWPACK to daily values and select nearest model grid cells."""
    station_x = xr.DataArray(
        snowpack_sd["x"].values,
        dims="station",
        coords={"station": snowpack_sd["station"].values},
    )
    station_y = xr.DataArray(
        snowpack_sd["y"].values,
        dims="station",
        coords={"station": snowpack_sd["station"].values},
    )

    snowpack_daily = snowpack_sd.resample(time="1D").mean()
    fsm_station = fsm_sd.sel(x=station_x, y=station_y, method="nearest")
    senorge_station = senorge_sd.sel(x=station_x, y=station_y, method="nearest")

    snowpack_daily, fsm_station, senorge_station = xr.align(
        snowpack_daily,
        fsm_station,
        senorge_station,
        join="inner",
    )
    return snowpack_daily, fsm_station, senorge_station


def nearest_grid_indices(grid_values: np.ndarray, point_values: np.ndarray) -> np.ndarray:
    """Return nearest 1D grid indices for each point value."""
    grid_values = np.asarray(grid_values)
    point_values = np.asarray(point_values)
    insertion = np.searchsorted(grid_values, point_values)
    insertion = np.clip(insertion, 1, len(grid_values) - 1)
    left = grid_values[insertion - 1]
    right = grid_values[insertion]
    return np.where(np.abs(point_values - left) <= np.abs(point_values - right), insertion - 1, insertion)


def extract_model_station_batch(
    nc_path: Path,
    variable_name: str,
    station_x: np.ndarray,
    station_y: np.ndarray,
    station_ids: np.ndarray,
    time_start: str,
    time_end: str,
    unit_divisor: float,
) -> xr.DataArray:
    """Extract daily model snow depth at station-nearest grid cells using netCDF4."""
    with Dataset(nc_path) as ds:
        grid_x = np.asarray(ds.variables["x"][:])
        grid_y = np.asarray(ds.variables["y"][:])
        x_idx = nearest_grid_indices(grid_x, station_x)
        y_idx = nearest_grid_indices(grid_y, station_y)

        time_var = ds.variables["time"]
        dates = np.array(
            num2date(
                time_var[:],
                time_var.units,
                getattr(time_var, "calendar", "standard"),
                only_use_cftime_datetimes=False,
                only_use_python_datetimes=True,
            ),
            dtype="datetime64[ns]",
        )
        time_mask = (dates >= np.datetime64(time_start)) & (dates <= np.datetime64(time_end))
        time_idx = np.where(time_mask)[0]
        selected_dates = dates[time_idx]

        values = np.empty((len(station_ids), len(time_idx)), dtype="float32")
        variable = ds.variables[variable_name]
        for row, (yi, xi) in enumerate(zip(y_idx, x_idx)):
            point_values = variable[time_idx, int(yi), int(xi)]
            values[row, :] = np.ma.filled(point_values, np.nan) / unit_divisor

    return xr.DataArray(
        values,
        dims=("station", "time"),
        coords={"station": station_ids, "time": selected_dates},
        name=variable_name,
        attrs={"units": "m"},
    )


def station_metrics(model: xr.DataArray, reference: xr.DataArray) -> xr.Dataset:
    """Compute bias, RMSE and Pearson correlation for each station."""
    model, reference = xr.align(model, reference, join="inner")
    valid = np.isfinite(model) & np.isfinite(reference)

    model_valid = model.where(valid)
    reference_valid = reference.where(valid)
    diff = model_valid - reference_valid

    n = valid.sum("time")
    bias = diff.mean("time", skipna=True)
    rmse = np.sqrt((diff**2).mean("time", skipna=True))

    model_anom = model_valid - model_valid.mean("time", skipna=True)
    reference_anom = reference_valid - reference_valid.mean("time", skipna=True)
    cov = (model_anom * reference_anom).sum("time", skipna=True)
    model_var = (model_anom**2).sum("time", skipna=True)
    reference_var = (reference_anom**2).sum("time", skipna=True)
    corr = cov / np.sqrt(model_var * reference_var)
    corr = corr.where(n >= 2)

    return xr.Dataset({"n": n, "bias": bias, "rmse": rmse, "corr": corr})


def compute_seasonal_metrics(
    snowpack_daily: xr.DataArray,
    fsm_station: xr.DataArray,
    senorge_station: xr.DataArray,
    seasons: list[Season] = DEFAULT_SEASONS,
) -> pd.DataFrame:
    """Compute station-wise seasonal metrics for FSM2 and seNorge."""
    rows = []

    for season in seasons:
        ref = snowpack_daily.sel(time=slice(season.start, season.end))
        fsm = fsm_station.sel(time=slice(season.start, season.end))
        senorge = senorge_station.sel(time=slice(season.start, season.end))

        fsm_metrics = station_metrics(fsm, ref)
        senorge_metrics = station_metrics(senorge, ref)

        df = pd.DataFrame(
            {
                "season": season.name,
                "season_label": season.label,
                "season_start": season.start,
                "season_end": season.end,
                "station": ref["station"].values,
                "elev": ref["elev"].values,
                "lon": ref["lon"].values,
                "lat": ref["lat"].values,
                "n_days": fsm_metrics["n"].values,
                "bias_fsm": fsm_metrics["bias"].values,
                "rmse_fsm": fsm_metrics["rmse"].values,
                "corr_fsm": fsm_metrics["corr"].values,
                "bias_senorge": senorge_metrics["bias"].values,
                "rmse_senorge": senorge_metrics["rmse"].values,
                "corr_senorge": senorge_metrics["corr"].values,
            }
        )
        df["rmse_diff_senorge_minus_fsm"] = df["rmse_senorge"] - df["rmse_fsm"]
        df["best_model_rmse"] = np.select(
            [
                df["rmse_fsm"] < df["rmse_senorge"],
                df["rmse_senorge"] < df["rmse_fsm"],
            ],
            ["FSM2", "seNorge"],
            default="Tie",
        )
        df["best_model_abs_bias"] = np.select(
            [
                df["bias_fsm"].abs() < df["bias_senorge"].abs(),
                df["bias_senorge"].abs() < df["bias_fsm"].abs(),
            ],
            ["FSM2", "seNorge"],
            default="Tie",
        )
        rows.append(df)

    return pd.concat(rows, ignore_index=True)


def summarize_metrics(metrics_df: pd.DataFrame) -> pd.DataFrame:
    """Create a compact season/model summary table."""
    records = []
    model_specs = [("FSM2", "fsm"), ("seNorge", "senorge")]

    for (season, season_label), season_df in metrics_df.groupby(["season", "season_label"], sort=False):
        for model_name, suffix in model_specs:
            records.append(
                {
                    "season": season,
                    "season_label": season_label,
                    "model": model_name,
                    "n_stations": int(season_df[f"rmse_{suffix}"].notna().sum()),
                    "mean_bias": season_df[f"bias_{suffix}"].mean(),
                    "median_bias": season_df[f"bias_{suffix}"].median(),
                    "mean_abs_bias": season_df[f"bias_{suffix}"].abs().mean(),
                    "mean_rmse": season_df[f"rmse_{suffix}"].mean(),
                    "median_rmse": season_df[f"rmse_{suffix}"].median(),
                    "mean_corr": season_df[f"corr_{suffix}"].mean(),
                    "median_corr": season_df[f"corr_{suffix}"].median(),
                }
            )

    return pd.DataFrame(records)


def compute_daily_error_evolution(
    snowpack_daily: xr.DataArray,
    fsm_station: xr.DataArray,
    senorge_station: xr.DataArray,
) -> pd.DataFrame:
    """Compute daily spatial mean bias and RMSE across all stations."""
    records = []
    for model_name, model in [("FSM2", fsm_station), ("seNorge", senorge_station)]:
        model, ref = xr.align(model, snowpack_daily, join="inner")
        diff = model - ref
        daily = xr.Dataset(
            {
                "mean_bias": diff.mean("station", skipna=True),
                "rmse": np.sqrt((diff**2).mean("station", skipna=True)),
                "mean_snow_depth": model.mean("station", skipna=True),
                "reference_mean_snow_depth": ref.mean("station", skipna=True),
            }
        ).to_dataframe()
        daily = daily.reset_index()
        daily["model"] = model_name
        records.append(daily)

    return pd.concat(records, ignore_index=True)


def daily_error_aggregate(
    model: xr.DataArray,
    reference: xr.DataArray,
    model_name: str,
) -> pd.DataFrame:
    """Aggregate daily error sums over stations for a station batch."""
    model, reference = xr.align(model, reference, join="inner")
    valid = np.isfinite(model) & np.isfinite(reference)
    diff = (model - reference).where(valid)

    ds = xr.Dataset(
        {
            "n": valid.sum("station"),
            "diff_sum": diff.sum("station", skipna=True),
            "diff2_sum": (diff**2).sum("station", skipna=True),
            "model_sum": model.where(valid).sum("station", skipna=True),
            "reference_sum": reference.where(valid).sum("station", skipna=True),
        }
    ).to_dataframe()
    ds = ds.reset_index()
    ds["model"] = model_name
    return ds


def finalize_daily_error_evolution(aggregate_df: pd.DataFrame) -> pd.DataFrame:
    """Convert batch-level daily sums to daily spatial mean bias and RMSE."""
    grouped = (
        aggregate_df.groupby(["time", "model"], as_index=False)[
            ["n", "diff_sum", "diff2_sum", "model_sum", "reference_sum"]
        ]
        .sum()
        .sort_values(["model", "time"])
    )
    grouped["mean_bias"] = grouped["diff_sum"] / grouped["n"]
    grouped["rmse"] = np.sqrt(grouped["diff2_sum"] / grouped["n"])
    grouped["mean_snow_depth"] = grouped["model_sum"] / grouped["n"]
    grouped["reference_mean_snow_depth"] = grouped["reference_sum"] / grouped["n"]
    return grouped[
        [
            "time",
            "model",
            "mean_bias",
            "rmse",
            "mean_snow_depth",
            "reference_mean_snow_depth",
            "n",
        ]
    ]


def compute_chunked_outputs(
    snowpack_sd: xr.DataArray,
    batch_size: int = 250,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute seasonal and daily outputs in station batches to keep memory low."""
    n_stations = snowpack_sd.sizes["station"]
    metric_frames = []
    daily_aggregate_frames = []

    for start in range(0, n_stations, batch_size):
        stop = min(start + batch_size, n_stations)
        print(f"Processing stations {start + 1}-{stop} of {n_stations}...", flush=True)

        snowpack_batch = snowpack_sd.isel(station=slice(start, stop))
        snowpack_daily = snowpack_batch.resample(time="1D").mean().load()

        fsm_station = extract_model_station_batch(
            PROCESSED / "fsm2_sd_all.nc",
            "snow_depth",
            snowpack_batch["x"].values,
            snowpack_batch["y"].values,
            snowpack_batch["station"].values,
            "2024-09-01",
            "2025-06-01",
            1000.0,
        )
        senorge_station = extract_model_station_batch(
            PROCESSED / "senorge_all.nc",
            "snow_depth",
            snowpack_batch["x"].values,
            snowpack_batch["y"].values,
            snowpack_batch["station"].values,
            "2024-09-01",
            "2025-06-01",
            100.0,
        )
        snowpack_daily, fsm_station, senorge_station = xr.align(
            snowpack_daily,
            fsm_station,
            senorge_station,
            join="inner",
        )

        metric_frames.append(
            compute_seasonal_metrics(
                snowpack_daily,
                fsm_station,
                senorge_station,
            )
        )
        daily_aggregate_frames.append(daily_error_aggregate(fsm_station, snowpack_daily, "FSM2"))
        daily_aggregate_frames.append(daily_error_aggregate(senorge_station, snowpack_daily, "seNorge"))

    metrics_df = pd.concat(metric_frames, ignore_index=True)
    daily_df = finalize_daily_error_evolution(pd.concat(daily_aggregate_frames, ignore_index=True))
    return metrics_df, daily_df


def plot_seasonal_boxplots(metrics_df: pd.DataFrame, output_path: Path) -> None:
    season_order = ["early", "midwinter", "melt", "full"]
    labels = ["Early", "Mid-winter", "Melt", "Full"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    metric_specs = [
        ("bias", "Bias [m]"),
        ("rmse", "RMSE [m]"),
        ("corr", "Correlation [-]"),
    ]

    for ax, (metric, ylabel) in zip(axes, metric_specs):
        positions = []
        box_data = []
        tick_positions = []

        for idx, season in enumerate(season_order):
            season_df = metrics_df[metrics_df["season"] == season]
            positions.extend([idx * 3 + 1, idx * 3 + 2])
            tick_positions.append(idx * 3 + 1.5)
            box_data.append(season_df[f"{metric}_fsm"].dropna().values)
            box_data.append(season_df[f"{metric}_senorge"].dropna().values)

        bp = ax.boxplot(
            box_data,
            positions=positions,
            widths=0.7,
            patch_artist=True,
            showfliers=False,
        )
        for patch, color in zip(bp["boxes"], ["#4c78a8", "#f58518"] * len(season_order)):
            patch.set_facecolor(color)
            patch.set_alpha(0.75)
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(labels, rotation=20)
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel.split(" [")[0])
        ax.grid(True, axis="y", alpha=0.25)

    axes[0].legend(
        [bp["boxes"][0], bp["boxes"][1]],
        ["FSM2", "seNorge"],
        loc="best",
        framealpha=0.9,
    )
    fig.suptitle("Seasonal Model Performance Against SNOWPACK", fontweight="bold")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_seasonal_rmse_maps(metrics_df: pd.DataFrame, output_path: Path) -> None:
    seasons = [
        ("early", "Early season"),
        ("midwinter", "Mid-winter"),
        ("melt", "Melt season"),
        ("full", "Full season"),
    ]
    vmax = metrics_df["rmse_diff_senorge_minus_fsm"].abs().quantile(0.98)
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex=True, sharey=True)

    for ax, (season, title) in zip(axes.ravel(), seasons):
        season_df = metrics_df[metrics_df["season"] == season]
        sc = ax.scatter(
            season_df["lon"],
            season_df["lat"],
            c=season_df["rmse_diff_senorge_minus_fsm"],
            s=10,
            cmap="RdBu",
            vmin=-vmax,
            vmax=vmax,
            linewidth=0,
            alpha=0.9,
            rasterized=True,
        )
        ax.set_title(title)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.grid(True, alpha=0.15)

    cbar = fig.colorbar(sc, ax=axes.ravel().tolist(), shrink=0.88, pad=0.02)
    cbar.set_label("RMSE seNorge - FSM2 [m]\nPositive values: FSM2 lower RMSE")
    fig.suptitle("Seasonal Spatial Difference in RMSE", fontweight="bold")
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_daily_error_evolution(daily_df: pd.DataFrame, output_path: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    colors = {"FSM2": "#4c78a8", "seNorge": "#f58518"}

    for model_name, model_df in daily_df.groupby("model"):
        axes[0].plot(model_df["time"], model_df["mean_bias"], label=model_name, color=colors[model_name])
        axes[1].plot(model_df["time"], model_df["rmse"], label=model_name, color=colors[model_name])

    for ax in axes:
        ax.axhline(0, color="black", linewidth=0.8, alpha=0.4)
        ax.axvspan(pd.Timestamp("2024-09-01"), pd.Timestamp("2024-12-31"), color="#6baed6", alpha=0.10)
        ax.axvspan(pd.Timestamp("2025-01-01"), pd.Timestamp("2025-03-31"), color="#9ecae1", alpha=0.13)
        ax.axvspan(pd.Timestamp("2025-04-01"), pd.Timestamp("2025-06-01"), color="#fdd0a2", alpha=0.20)
        ax.grid(True, alpha=0.25)
        ax.legend(loc="best")

    axes[0].set_ylabel("Mean bias [m]")
    axes[0].set_title("Daily Spatial Mean Bias")
    axes[1].set_ylabel("RMSE [m]")
    axes[1].set_title("Daily Spatial RMSE")
    axes[1].set_xlabel("Time")
    fig.suptitle("Seasonal Evolution of Model Errors", fontweight="bold")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    print("Loading SNOWPACK data...", flush=True)
    snowpack_sd = load_snowpack_depth()

    print("Computing seasonal station metrics and daily error evolution...", flush=True)
    metrics_df, daily_df = compute_chunked_outputs(snowpack_sd)
    summary_df = summarize_metrics(metrics_df)

    metrics_path = PROCESSED / "seasonal_model_metrics.csv"
    summary_path = PROCESSED / "seasonal_model_summary.csv"
    metrics_df.to_csv(metrics_path, index=False)
    summary_df.to_csv(summary_path, index=False)

    daily_path = PROCESSED / "seasonal_daily_error_evolution.csv"
    daily_df.to_csv(daily_path, index=False)

    print("Creating figures...")
    plot_seasonal_boxplots(metrics_df, PROCESSED / "Figure_09_Seasonal_Model_Performance.png")
    plot_seasonal_rmse_maps(metrics_df, PROCESSED / "Figure_10_Seasonal_RMSE_Difference_Maps.png")
    plot_daily_error_evolution(daily_df, PROCESSED / "Figure_11_Seasonal_Error_Evolution.png")

    print("\nSaved outputs:")
    for path in [
        metrics_path,
        summary_path,
        daily_path,
        PROCESSED / "Figure_09_Seasonal_Model_Performance.png",
        PROCESSED / "Figure_10_Seasonal_RMSE_Difference_Maps.png",
        PROCESSED / "Figure_11_Seasonal_Error_Evolution.png",
    ]:
        print(f"- {path.relative_to(ROOT)}")

    print("\nSeasonal summary:")
    print(summary_df.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
