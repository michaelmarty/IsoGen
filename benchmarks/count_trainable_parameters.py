"""Print trainable parameter counts for the neural models used by the timing scripts.

Run from the repository root with:
    python benchmarks/count_trainable_parameters.py
"""

from pathlib import Path

import torch


MODEL_NAMES = (
    "isogenpep_model_16",
    "isogenpep_model_64",
    "isogenpep_model_128",
    "isogenmass_model_8",
    "isogenmass_model_32",
    "isogenmass_model_64",
    "isogenmass_model_128",
    "isogenmass_model_1024",
    "isogenrna_model_64",
    "isogenrna_model_128",
    "isogen_rnaveragine_model32",
    "isogen_rnaveragine_model64",
    "isogen_rnaveragine_model128",
)


def count_parameters(model_path: Path) -> int:
    weights = torch.load(model_path, map_location="cpu", weights_only=True)
    return sum(tensor.numel() for tensor in weights.values())


def main() -> None:
    model_directory = Path(__file__).resolve().parents[1] / "isogen" / "models"
    for name in MODEL_NAMES:
        parameter_count = count_parameters(model_directory / f"{name}.pth")
        print(f"{name}: {parameter_count:,} trainable parameters")

    print("FFT and BRAIN are native algorithms and have no trainable parameters.")


if __name__ == "__main__":
    main()
