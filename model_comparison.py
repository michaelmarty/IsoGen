"""Compare IsoGen mass-model depth and input encodings on held-out masses.

Edit the experiment parameters below, then run ``python model_comparison.py``.
The FFT peptide averagine calculation supplies the reference distributions.
"""

from __future__ import annotations

import csv
import math
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset

from isogen.isogen_tools import mass_to_vector
from isogen.isogenwrapper import fft_gen_isodist


# Experiment parameters (edit these values as needed).
NUM_TRAINING_MASSES = 10_000
NUM_VALIDATION_MASSES = 2_000
NUM_TEST_MASSES = 10_000
MIN_MASS_DA = 100.0
MAX_MASS_DA = 100_000.0
ISOTOPE_VECTOR_LENGTH = 128
HIDDEN_LAYER_WIDTH = 128
HIDDEN_LAYER_DEPTHS = (0, 1, 2, 3)
ENCODING_STRATEGIES = ("Current", "Normalized", "Log", "Direct")
STANDARD_MODEL_DEPTH = 2
EPOCHS = 10
BATCH_SIZE = 32
LEARNING_RATE = 0.001
RANDOM_SEED = 2026
RESULTS_CSV = Path("C:\\Users\\mm96978\\The University of Texas at Austin\\MartyLab - General\\Papers\\IsoGen\\Figures\\model_comparison_results.csv")
LOSS_CURVE_PNG = Path("C:\\Users\\mm96978\\The University of Texas at Austin\\MartyLab - General\\Papers\\IsoGen\\Figures\\model_comparison_loss_curves.png")
LOSS_CURVE_PDF = Path("C:\\Users\\mm96978\\The University of Texas at Austin\\MartyLab - General\\Papers\\IsoGen\\Figures\\model_comparison_loss_curves.pdf")


def encode_masses(masses: np.ndarray, strategy: str) -> np.ndarray:
    """Return one of the four reviewer-requested mass representations."""
    masses = np.asarray(masses, dtype=np.float32)
    if strategy == "Current":
        return np.asarray([mass_to_vector(float(mass)) for mass in masses], dtype=np.float32)
    if strategy == "Normalized":
        return ((masses - MIN_MASS_DA) / (MAX_MASS_DA - MIN_MASS_DA))[:, None]
    if strategy == "Log":
        return np.log10(masses)[:, None]
    if strategy == "Direct":
        return masses[:, None]
    raise ValueError(f"Unknown encoding strategy: {strategy}")


class MassModel(nn.Module):
    """IsoGen's standard Softsign/Sigmoid model with configurable depth."""

    def __init__(self, input_width: int, hidden_layers: int):
        super().__init__()
        layers: list[nn.Module] = []
        width = input_width
        for _ in range(hidden_layers):
            layers.extend((nn.Linear(width, HIDDEN_LAYER_WIDTH), nn.Softsign()))
            width = HIDDEN_LAYER_WIDTH
        layers.extend((nn.Linear(width, ISOTOPE_VECTOR_LENGTH), nn.Sigmoid()))
        self.network = nn.Sequential(*layers)

    def forward(self, inputs: torch.Tensor) -> torch.Tensor:
        return self.network(inputs)


@dataclass
class Result:
    comparison: str
    model: str
    hidden_layers: int
    input_width: int
    parameters: int
    train_seconds: float
    mse: float
    rmse: float
    mae: float
    mean_vector_sse: float
    mean_cosine_similarity: float


@dataclass
class History:
    training_loss: list[float]
    validation_loss: list[float]
    training_css: list[float]
    validation_css: list[float]


def generate_dataset(rng: np.random.Generator, size: int, label: str) -> tuple[np.ndarray, np.ndarray]:
    masses = rng.uniform(MIN_MASS_DA, MAX_MASS_DA, size).astype(np.float32)
    print(f"Generating {label} FFT targets for {size:,} masses...")
    targets = np.asarray(
        [fft_gen_isodist(float(mass), type="PEPTIDE", isolen=ISOTOPE_VECTOR_LENGTH) for mass in masses],
        dtype=np.float32,
    )
    return masses, targets


def train_and_evaluate(
    train_masses: np.ndarray,
    train_targets: np.ndarray,
    validation_masses: np.ndarray,
    validation_targets: np.ndarray,
    test_masses: np.ndarray,
    test_targets: np.ndarray,
    encoding: str,
    hidden_layers: int,
    device: torch.device,
) -> tuple[Result, History]:
    train_inputs = torch.from_numpy(encode_masses(train_masses, encoding))
    validation_inputs = torch.from_numpy(encode_masses(validation_masses, encoding))
    test_inputs = torch.from_numpy(encode_masses(test_masses, encoding))
    train_outputs = torch.from_numpy(train_targets)
    validation_outputs = torch.from_numpy(validation_targets)
    test_outputs = torch.from_numpy(test_targets)

    torch.manual_seed(RANDOM_SEED)
    if device.type == "cuda":
        torch.cuda.manual_seed_all(RANDOM_SEED)
    model = MassModel(train_inputs.shape[1], hidden_layers).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=LEARNING_RATE)
    loss_function = nn.MSELoss()
    generator = torch.Generator().manual_seed(RANDOM_SEED)
    loader = DataLoader(
        TensorDataset(train_inputs, train_outputs),
        batch_size=BATCH_SIZE,
        shuffle=True,
        generator=generator,
        pin_memory=device.type == "cuda",
    )

    train_seconds = 0.0
    history = History([], [], [], [])
    for epoch in range(EPOCHS):
        started = time.perf_counter()
        model.train()
        total_loss = 0.0
        for inputs, targets in loader:
            inputs = inputs.to(device, non_blocking=True)
            targets = targets.to(device, non_blocking=True)
            optimizer.zero_grad(set_to_none=True)
            loss = loss_function(model(inputs), targets)
            loss.backward()
            optimizer.step()
            total_loss += loss.item() * len(inputs)
        train_seconds += time.perf_counter() - started
        training_mse = total_loss / len(loader.dataset)
        _, training_css = calculate_mse_and_css(model, train_inputs, train_outputs, device)
        validation_mse, validation_css = calculate_mse_and_css(
            model, validation_inputs, validation_outputs, device
        )
        history.training_loss.append(training_mse)
        history.validation_loss.append(validation_mse)
        history.training_css.append(training_css)
        history.validation_css.append(validation_css)
        print(
            f"  encoding={encoding:<10} depth={hidden_layers} "
            f"epoch={epoch + 1:>2}/{EPOCHS} train_mse={training_mse:.7g} "
            f"validation_mse={validation_mse:.7g} validation_css={validation_css:.5f}"
        )

    model.eval()
    squared_error_sum = 0.0
    absolute_error_sum = 0.0
    vector_sse_sum = 0.0
    cosine_similarity_sum = 0.0
    value_count = 0
    with torch.no_grad():
        for start in range(0, len(test_inputs), BATCH_SIZE):
            inputs = test_inputs[start : start + BATCH_SIZE].to(device)
            targets = test_outputs[start : start + BATCH_SIZE].to(device)
            predictions = model(inputs)
            errors = predictions - targets
            squared = errors.square()
            squared_error_sum += squared.sum().item()
            absolute_error_sum += errors.abs().sum().item()
            vector_sse_sum += squared.sum(dim=1).sum().item()
            cosine_similarity_sum += nn.functional.cosine_similarity(
                predictions, targets, dim=1
            ).sum().item()
            value_count += errors.numel()

    mse = squared_error_sum / value_count
    return (
        Result(
            comparison="",
            model=encoding,
            hidden_layers=hidden_layers,
            input_width=train_inputs.shape[1],
            parameters=sum(parameter.numel() for parameter in model.parameters()),
            train_seconds=train_seconds,
            mse=mse,
            rmse=math.sqrt(mse),
            mae=absolute_error_sum / value_count,
            mean_vector_sse=vector_sse_sum / len(test_inputs),
            mean_cosine_similarity=cosine_similarity_sum / len(test_inputs),
        ),
        history,
    )


def calculate_mse_and_css(
    model: nn.Module,
    inputs: torch.Tensor,
    targets: torch.Tensor,
    device: torch.device,
) -> tuple[float, float]:
    model.eval()
    squared_error_sum = 0.0
    cosine_similarity_sum = 0.0
    with torch.no_grad():
        for start in range(0, len(inputs), BATCH_SIZE):
            batch_inputs = inputs[start : start + BATCH_SIZE].to(device)
            batch_targets = targets[start : start + BATCH_SIZE].to(device)
            predictions = model(batch_inputs)
            squared_error_sum += (predictions - batch_targets).square().sum().item()
            cosine_similarity_sum += nn.functional.cosine_similarity(
                predictions, batch_targets, dim=1
            ).sum().item()
    return squared_error_sum / targets.numel(), cosine_similarity_sum / len(targets)


def print_results(title: str, results: list[Result]) -> None:
    print(f"\n{title}")
    print("model        depth inputs   parameters  train (s)       MSE      RMSE       MAE  mean SSE  cosine")
    for result in results:
        print(
            f"{result.model:<12} {result.hidden_layers:>5} {result.input_width:>6} "
            f"{result.parameters:>12,} {result.train_seconds:>8.2f} "
            f"{result.mse:>9.3g} {result.rmse:>9.3g} {result.mae:>9.3g} "
            f"{result.mean_vector_sse:>9.3g} {result.mean_cosine_similarity:>7.4f}"
        )


def plot_loss_curves(
    histories: dict[tuple[int, str], tuple[Result, History]],
) -> None:
    figure, axes = plt.subplots(2, 1, figsize=(11, 10), sharex=True)
    epochs = range(1, EPOCHS + 1)
    for (depth, encoding), (_, history) in histories.items():
        label = f"{depth} Hidden/Vector" if encoding == "Current" else f"{depth} Hidden/{encoding}"
        training_line = axes[0].plot(epochs, history.training_loss, label=f"{label} train")[0]
        color = training_line.get_color()
        axes[0].plot(
            epochs,
            history.validation_loss,
            linestyle="--",
            color=color,
            label=f"{label} validation",
        )
        axes[1].plot(epochs, history.training_css, color=color)
        axes[1].plot(epochs, history.validation_css, linestyle="--", color=color)
    axes[0].set(ylabel="Mean Squared Error", title="Training and Validation Loss")
    axes[0].set_yscale("log")
    axes[0].legend(ncol=2, fontsize=8)
    axes[1].set(
        xlabel="Epoch",
        ylabel="Mean Cosine Similarity Score",
        title="Training and Validation CSS",
        ylim=(0.0, 1.01),
    )
    axes[0].text(0.98, 0.96, "A", transform=axes[0].transAxes, ha="right", va="top", fontsize=16, fontweight="bold")
    axes[1].text(0.98, 0.04, "B", transform=axes[1].transAxes, ha="right", va="bottom", fontsize=16, fontweight="bold")
    for axis in axes:
        axis.grid(True, which="both", alpha=0.25)
    # figure.suptitle("IsoGen Model Comparison")
    figure.tight_layout()
    figure.savefig(LOSS_CURVE_PNG, dpi=300, bbox_inches="tight")
    figure.savefig(LOSS_CURVE_PDF, bbox_inches="tight")
    plt.close(figure)


def validate_parameters() -> None:
    if min(NUM_TRAINING_MASSES, NUM_VALIDATION_MASSES, NUM_TEST_MASSES) < 1:
        raise ValueError("Training, validation, and test sample counts must be positive")
    if not 0 < MIN_MASS_DA < MAX_MASS_DA:
        raise ValueError("Mass limits must satisfy 0 < MIN_MASS_DA < MAX_MASS_DA")
    if ISOTOPE_VECTOR_LENGTH < 1 or HIDDEN_LAYER_WIDTH < 1:
        raise ValueError("Model widths must be positive")
    if EPOCHS < 1 or BATCH_SIZE < 1:
        raise ValueError("EPOCHS and BATCH_SIZE must be positive")
    if any(depth not in (0, 1, 2, 3) for depth in HIDDEN_LAYER_DEPTHS):
        raise ValueError("HIDDEN_LAYER_DEPTHS may contain only 0, 1, 2, or 3")
    unknown = set(ENCODING_STRATEGIES) - {"Current", "Normalized", "Log", "Direct"}
    if unknown:
        raise ValueError(f"Unknown encodings: {sorted(unknown)}")


def main() -> None:
    validate_parameters()
    rng = np.random.default_rng(RANDOM_SEED)
    train_masses, train_targets = generate_dataset(rng, NUM_TRAINING_MASSES, "training")
    validation_masses, validation_targets = generate_dataset(
        rng, NUM_VALIDATION_MASSES, "validation"
    )
    test_masses, test_targets = generate_dataset(rng, NUM_TEST_MASSES, "held-out test")
    device = torch.device(
        "cuda" if torch.cuda.is_available() else "mps" if torch.backends.mps.is_available() else "cpu"
    )
    print(f"Training on {device}; all models use Adam, MSE, Softsign hidden activations, and Sigmoid output.")

    cache: dict[tuple[int, str], tuple[Result, History]] = {}

    def run(depth: int, encoding: str) -> tuple[Result, History]:
        key = (depth, encoding)
        if key not in cache:
            cache[key] = train_and_evaluate(
                train_masses,
                train_targets,
                validation_masses,
                validation_targets,
                test_masses,
                test_targets,
                encoding,
                depth,
                device,
            )
        return cache[key]

    depth_results = []
    for depth in HIDDEN_LAYER_DEPTHS:
        result, _ = run(depth, "Current")
        depth_results.append(Result(**{**result.__dict__, "comparison": "depth", "model": f"{depth} hidden"}))

    encoding_results = []
    for encoding in ENCODING_STRATEGIES:
        result, _ = run(STANDARD_MODEL_DEPTH, encoding)
        encoding_results.append(Result(**{**result.__dict__, "comparison": "encoding"}))

    print_results("Depth comparison (current five-value mass encoding)", depth_results)
    print_results(f"Encoding comparison ({STANDARD_MODEL_DEPTH} hidden layers)", encoding_results)
    plot_loss_curves(cache)

    all_results = depth_results + encoding_results
    with RESULTS_CSV.open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=Result.__dataclass_fields__)
        writer.writeheader()
        writer.writerows(result.__dict__ for result in all_results)
    print(f"\nSaved results to {RESULTS_CSV.resolve()}")
    print(f"Saved loss curves to {LOSS_CURVE_PNG.resolve()} and {LOSS_CURVE_PDF.resolve()}")


if __name__ == "__main__":
    main()
