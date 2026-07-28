"""
main.py

Calcula:
    - Volumen del material
    - Volumen de la concavidad (capacidad)
    - Visualización 3D de voxeles

"""

from time import perf_counter

from config import MODEL_FILE, POINTS_PER_BATCH, VOXEL_SIZE_MM
from voxelizer import VolumeCalculator


def main():

    inicio = perf_counter()

    calc = VolumeCalculator(
        MODEL_FILE,
        VOXEL_SIZE_MM,
        POINTS_PER_BATCH
    )

    print("\n==============================")
    print("   ANALISIS DEL MODELO 3D")
    print("==============================\n")

    # -------------------------
    # MATERIAL
    # -------------------------

    material_volume, voxel_count = calc.calculate()

    # -------------------------
    # CONCAVIDAD
    # -------------------------

    print("\nBuscando concavidades...")

    try:

        cavity_volume = calc.calculate_void()

    except AttributeError:

        cavity_volume = None

        print("calculate_void() todavía no implementado.")

    tiempo = perf_counter() - inicio

    print("\n===================================")
    print("RESULTADOS")
    print("===================================")

    print(f"Voxel utilizado      : {VOXEL_SIZE_MM} mm")

    print(f"Voxeles internos     : {voxel_count:,}")

    print(f"\nVolumen material")

    print(f"    {material_volume:.2f} mm³")

    print(f"    {material_volume/1000:.2f} cm³")

    if cavity_volume is not None:

        print(f"\nVolumen concavidad")

        print(f"    {cavity_volume:.2f} mm³")

        print(f"    {cavity_volume/1000:.2f} cm³")

        print(f"\nVolumen total")

        print(f"    {(material_volume+cavity_volume)/1000:.2f} cm³")

    print(f"\nTiempo: {tiempo:.2f} s")

    print("===================================")

    # -------------------------
    # VISUALIZACION
    # -------------------------

    try:

        calc.show_voxels()

    except AttributeError:

        print("\nshow_voxels() todavía no implementado.")


if __name__ == "__main__":
    main()