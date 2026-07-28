"""Cálculo aproximado de volumen por voxelización para mallas 3D."""

from pathlib import Path
from time import perf_counter

import numpy as np
import trimesh


class VolumeCalculator:
    """Cuenta voxeles cuyos centros se encuentran dentro de una malla cerrada."""

    def __init__(
        self, filename: str, voxel_size_mm: float, points_per_batch: int = 20_000
    ):
        self.filename = Path(filename)
        self.voxel_size_mm = float(voxel_size_mm)
        self.points_per_batch = int(points_per_batch)

        if self.voxel_size_mm <= 0:
            raise ValueError("El tamaño del voxel debe ser mayor que cero.")
        if self.points_per_batch <= 0:
            raise ValueError("POINTS_PER_BATCH must be greater than zero.")
        if not self.filename.is_file():
            raise FileNotFoundError(
                f"No se encontró '{self.filename}'. Coloca el archivo OBJ junto a main.py "
                "o cambia MODEL_FILE en config.py."
            )

        # Un archivo MTL describe colores/materiales, pero no contiene la forma 3D.
        # A veces se guarda por error con extensión .obj; detectarlo evita confusión.
        with self.filename.open("r", encoding="utf-8", errors="ignore") as model_file:
            first_part = model_file.read(300)
        if "newmtl " in first_part and not any(
            line.startswith("v ") or line.startswith("f ")
            for line in first_part.splitlines()
        ):
            raise ValueError(
                f"'{self.filename.name}' contiene materiales (MTL), no la malla OBJ. "
                "Exporta nuevamente el cuerpo como OBJ o STL."
            )

        loaded = trimesh.load(self.filename, force="mesh")
        if isinstance(loaded, trimesh.Scene):
            loaded = trimesh.util.concatenate(tuple(loaded.geometry.values()))
        self.mesh = loaded

        if self.mesh.is_empty:
            raise ValueError("El archivo no contiene una malla con caras para calcular.")

    @staticmethod
    def _show_progress(done: int, total: int, start_time: float) -> None:
        """Show a progress bar that refreshes in the terminal."""
        ratio = min(done / total, 1.0) if total else 1.0
        width = 30
        filled = int(ratio * width)
        bar = "#" * filled + "-" * (width - filled)
        elapsed = perf_counter() - start_time

        if done and elapsed > 0:
            rate = done / elapsed
            remaining = (total - done) / rate
            status = f"{ratio * 100:6.2f}% | {done:,}/{total:,} voxels | ETA {remaining:,.0f}s"
        else:
            status = "  0.00% | preparing calculation"

        print(f"\r[{bar}] {status}", end="", flush=True)
        if done >= total:
            print()

    def calculate(self) -> tuple[float, int]:
        """Devuelve (volumen_mm3, cantidad_de_voxeles_dentro)."""
        if not self.mesh.is_watertight:
            print("ADVERTENCIA: la malla no está cerrada (watertight).")
            print("El resultado puede no representar el volumen real.\n")

        minimum, maximum = self.mesh.bounds
        x_values = np.arange(
            minimum[0] + self.voxel_size_mm / 2,
            maximum[0],
            self.voxel_size_mm,
        )
        y_values = np.arange(
            minimum[1] + self.voxel_size_mm / 2,
            maximum[1],
            self.voxel_size_mm,
        )
        z_values = np.arange(
            minimum[2] + self.voxel_size_mm / 2,
            maximum[2],
            self.voxel_size_mm,
        )

        print("Caja envolvente (mm):")
        print(self.mesh.bounds)
        print(f"Tamaño de voxel: {self.voxel_size_mm} mm")

        # Process each layer in small batches. This keeps RAM use low and lets
        # the terminal update continuously instead of waiting for a full layer.
        grid_x, grid_y = np.meshgrid(x_values, y_values, indexing="ij")
        x_flat = grid_x.ravel()
        y_flat = grid_y.ravel()
        voxels_per_layer = x_flat.size
        total_voxels = voxels_per_layer * len(z_values)
        print(f"Voxeles por evaluar: {total_voxels:,}")
        print("Progreso:")

        total_inside = 0
        processed = 0
        start_time = perf_counter()
        self._show_progress(processed, total_voxels, start_time)

        for z_value in z_values:
            for start in range(0, voxels_per_layer, self.points_per_batch):
                end = min(start + self.points_per_batch, voxels_per_layer)
                points = np.column_stack(
                    (
                        x_flat[start:end],
                        y_flat[start:end],
                        np.full(end - start, z_value),
                    )
                )
                inside = self.mesh.contains(points)
                total_inside += int(np.count_nonzero(inside))
                processed += end - start
                self._show_progress(processed, total_voxels, start_time)

        volume_mm3 = total_inside * self.voxel_size_mm**3
        return volume_mm3, total_inside


# ============================================================
# NUEVAS FUNCIONES (EXPERIMENTALES)
# ============================================================

    def calculate_void(self):
        """
        Aproximación del volumen de la concavidad abierta hacia +Z.

        Para cada columna X,Y busca la primera intersección con el sólido.
        El espacio entre Zmax y esa primera superficie se considera capacidad.
        """

        minimum, maximum = self.mesh.bounds

        x_values = np.arange(
            minimum[0] + self.voxel_size_mm / 2,
            maximum[0],
            self.voxel_size_mm,
        )

        y_values = np.arange(
            minimum[1] + self.voxel_size_mm / 2,
            maximum[1],
            self.voxel_size_mm,
        )

        z_values = np.arange(
            maximum[2] - self.voxel_size_mm / 2,
            minimum[2],
            -self.voxel_size_mm,
        )

        cavity = 0

        for x in x_values:
            for y in y_values:

                started = False

                for z in z_values:

                    p = np.array([[x, y, z]])

                    inside = self.mesh.contains(p)[0]

                    if inside:
                        started = True

                    elif started:
                        cavity += 1

                    elif not started:
                        pass

        return cavity * self.voxel_size_mm**3


    def show_voxels(self):
        """
        Visualización sencilla de los centros de voxeles del material.
        Requiere:
            pip install pyvista
        """

        try:
            import pyvista as pv
        except ImportError:
            print("Instala pyvista con:")
            print("python -m pip install pyvista")
            return

        minimum, maximum = self.mesh.bounds

        x = np.arange(
            minimum[0] + self.voxel_size_mm / 2,
            maximum[0],
            self.voxel_size_mm,
        )

        y = np.arange(
            minimum[1] + self.voxel_size_mm / 2,
            maximum[1],
            self.voxel_size_mm,
        )

        z = np.arange(
            minimum[2] + self.voxel_size_mm / 2,
            maximum[2],
            self.voxel_size_mm,
        )

        pts = []

        for zz in z:
            X, Y = np.meshgrid(x, y, indexing="ij")

            points = np.column_stack(
                (
                    X.ravel(),
                    Y.ravel(),
                    np.full(X.size, zz),
                )
            )

            inside = self.mesh.contains(points)

            pts.extend(points[inside])

        cloud = pv.PolyData(np.array(pts))

        plotter = pv.Plotter()

        plotter.add_mesh(
            self.mesh,
            opacity=0.20,
            color="lightgray",
        )

        plotter.add_points(
            cloud,
            color="blue",
            point_size=6,
            render_points_as_spheres=True,
        )

        plotter.show()
