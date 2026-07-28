"""Example for measuring an irregular wound region.

Edit wound_boundary.csv first, then run: python main_wound_metrics.py
"""

import json
from pathlib import Path

from config import MODEL_FILE
from wound_metrics import calculate_wound_metrics


BOUNDARY_FILE = "wound_boundary.csv"
RESULT_FILE = Path("outputs") / "wound_metrics.json"


def show_progress(percent: int, label: str) -> None:
    """Display a simple loading bar for the analysis stages."""
    width = 30
    filled = int(percent / 100 * width)
    bar = "#" * filled + "-" * (width - filled)
    print(f"\r[{bar}] {percent:3d}%  {label:<42}", end="", flush=True)
    if percent >= 100:
        print()


def main() -> None:
    print("Starting wound analysis...")
    metrics = calculate_wound_metrics(MODEL_FILE, BOUNDARY_FILE, show_progress)
    show_progress(100, "Saving JSON report")
    RESULT_FILE.parent.mkdir(exist_ok=True)
    RESULT_FILE.write_text(json.dumps(metrics, indent=2), encoding="utf-8")

    print("\nWOUND METRICS")
    print("-------------")
    for name, value in metrics.items():
        print(f"{name}: {value}")
    print(f"\nSaved JSON report: {RESULT_FILE}")


if __name__ == "__main__":
    main()
