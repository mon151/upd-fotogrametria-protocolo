"""Run this file to calculate the approximate volume of the selected model.

Steps:
1. Put the OBJ or STL file next to this file.
2. Set its name in config.py.
3. Run: python main.py

Required packages (installed once):
    python -m pip install numpy trimesh rtree
"""

from time import perf_counter

from config import MODEL_FILE, POINTS_PER_BATCH, VOXEL_SIZE_MM
from voxelizer import VolumeCalculator


def main() -> None:
    start_time = perf_counter()
    calculator = VolumeCalculator(MODEL_FILE, VOXEL_SIZE_MM, POINTS_PER_BATCH)
    volume_mm3, voxel_count = calculator.calculate()
    elapsed_seconds = perf_counter() - start_time

    print("\n-------------------------")
    print("RESULT")
    print("-------------------------")
    print(f"Voxels inside the object: {voxel_count:,}")
    print(f"Approximate volume: {volume_mm3:.3f} mm3")
    print(f"Approximate volume: {volume_mm3 / 1000:.3f} cm3")
    print(f"Calculation time: {elapsed_seconds:.1f} seconds")
    print("-------------------------")


if __name__ == "__main__":
    main()
