"""Geometry-only measurements for an irregular wound region on a 3D mesh."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Callable

import numpy as np
import trimesh


def load_boundary_csv(filename: str | Path) -> np.ndarray:
    """Read ordered x_mm,y_mm,z_mm points that trace the wound edge."""
    rows: list[list[float]] = []
    with Path(filename).open(newline="", encoding="utf-8") as csv_file:
        for row in csv.DictReader(csv_file):
            rows.append([float(row["x_mm"]), float(row["y_mm"]), float(row["z_mm"])])
    boundary = np.asarray(rows, dtype=float)
    if len(boundary) < 3:
        raise ValueError(
            "wound_boundary.csv is empty. Add at least 3 ordered edge points "
            "with x_mm,y_mm,z_mm coordinates from the same 3D model."
        )
    return boundary


def _local_coordinates(points: np.ndarray, boundary: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Project 3D points to local u, v, height coordinates defined by the rim."""
    center = boundary.mean(axis=0)
    _, _, vectors = np.linalg.svd(boundary - center, full_matrices=False)
    axes = vectors[:2]
    normal = np.cross(axes[0], axes[1])
    normal /= np.linalg.norm(normal)
    relative = points - center
    local = np.column_stack((relative @ axes[0], relative @ axes[1], relative @ normal))
    return local, center, normal


def _inside_polygon(points: np.ndarray, polygon: np.ndarray) -> np.ndarray:
    """Return whether each 2D point lies inside an ordered, possibly irregular rim."""
    x, y = points[:, 0], points[:, 1]
    px, py = polygon[:, 0], polygon[:, 1]
    inside = np.zeros(len(points), dtype=bool)
    previous = len(polygon) - 1
    for current in range(len(polygon)):
        crosses = (py[current] > y) != (py[previous] > y)
        intersection_x = (px[previous] - px[current]) * (y - py[current]) / (
            py[previous] - py[current] + 1e-15
        ) + px[current]
        inside ^= crosses & (x < intersection_x)
        previous = current
    return inside


def _thin_plate_reference(rim_uv: np.ndarray, rim_height: np.ndarray, target_uv: np.ndarray) -> np.ndarray:
    """Smooth reference surface through an irregular 3D rim; no sphere is used."""
    distances = np.linalg.norm(rim_uv[:, None] - rim_uv[None, :], axis=2)
    kernel = distances**2 * np.log(distances + 1e-12)
    count = len(rim_uv)
    polynomial = np.column_stack((np.ones(count), rim_uv))
    system = np.block(
        [[kernel + np.eye(count) * 1e-10, polynomial], [polynomial.T, np.zeros((3, 3))]]
    )
    weights = np.linalg.solve(system, np.concatenate((rim_height, np.zeros(3))))

    target_distances = np.linalg.norm(target_uv[:, None] - rim_uv[None, :], axis=2)
    target_kernel = target_distances**2 * np.log(target_distances + 1e-12)
    return target_kernel @ weights[:count] + np.column_stack((np.ones(len(target_uv)), target_uv)) @ weights[count:]


def _polygon_area(polygon: np.ndarray) -> float:
    shifted = np.roll(polygon, -1, axis=0)
    return float(abs(np.sum(polygon[:, 0] * shifted[:, 1] - shifted[:, 0] * polygon[:, 1])) / 2)


def _feret_diameter(points: np.ndarray) -> float:
    distances = np.linalg.norm(points[:, None] - points[None, :], axis=2)
    return float(distances.max())


def calculate_wound_metrics(
    mesh_file: str | Path,
    boundary_file: str | Path,
    progress: Callable[[int, str], None] | None = None,
) -> dict:
    """Calculate generalized wound geometry from a mesh and manually traced rim."""
    report = progress or (lambda _percent, _label: None)
    report(10, "Loading 3D mesh")
    loaded = trimesh.load(mesh_file, force="mesh")
    mesh = trimesh.util.concatenate(tuple(loaded.geometry.values())) if isinstance(loaded, trimesh.Scene) else loaded
    if mesh.is_empty:
        raise ValueError("The mesh has no faces.")

    report(25, "Loading wound boundary")
    boundary = load_boundary_csv(boundary_file)
    report(40, "Building local reference coordinates")
    vertex_local, _, normal = _local_coordinates(mesh.vertices, boundary)
    rim_local, _, _ = _local_coordinates(boundary, boundary)
    rim_uv = rim_local[:, :2]
    vertex_inside = _inside_polygon(vertex_local[:, :2], rim_uv)
    if vertex_inside.sum() < 3:
        raise ValueError("No mesh vertices were found inside the supplied wound boundary.")

    report(55, "Estimating smooth reference surface")
    reference_height = _thin_plate_reference(rim_uv, rim_local[:, 2], vertex_local[:, :2])
    signed_depth = reference_height - vertex_local[:, 2]
    positive_total = np.maximum(signed_depth[vertex_inside], 0).sum()
    negative_total = np.maximum(-signed_depth[vertex_inside], 0).sum()
    orientation = 1.0 if positive_total >= negative_total else -1.0
    depth = np.maximum(orientation * signed_depth, 0)

    report(70, "Calculating depth, area and volume")
    faces = mesh.faces
    face_centers = mesh.triangles_center
    center_local, _, _ = _local_coordinates(face_centers, boundary)
    face_inside = _inside_polygon(center_local[:, :2], rim_uv)
    center_reference = _thin_plate_reference(rim_uv, rim_local[:, 2], center_local[:, :2])
    face_depth = np.maximum(orientation * (center_reference - center_local[:, 2]), 0)

    triangles = mesh.triangles
    projected_face_area = 0.5 * np.abs(np.einsum("ij,j->i", np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]), normal))
    selected_faces = face_inside & (face_depth > 0)
    selected_depth = depth[vertex_inside & (depth > 0)]
    if len(selected_depth) == 0:
        raise ValueError("The selected region has no measurable concavity relative to its rim.")

    rim_closed = np.vstack((boundary, boundary[0]))
    perimeter = float(np.linalg.norm(np.diff(rim_closed, axis=0), axis=1).sum())
    face_normals = mesh.face_normals[selected_faces]
    slopes = np.degrees(np.arccos(np.clip(np.abs(face_normals @ normal), 0, 1)))

    report(85, "Calculating wall and width profile")
    depth_profile: dict[str, float | None] = {}
    max_depth = float(selected_depth.max())
    for percentage in (25, 50, 75):
        target = max_depth * percentage / 100
        band = max(max_depth * 0.08, 0.25)
        values = vertex_local[vertex_inside & (np.abs(depth - target) <= band), :2]
        depth_profile[f"width_at_{percentage}_percent_depth_mm"] = (
            float(_feret_diameter(values)) if len(values) >= 2 else None
        )

    report(100, "Measurements complete")
    return {
        "units": "mm, mm2, mm3, degrees",
        "mesh_is_watertight": bool(mesh.is_watertight),
        "boundary_points": int(len(boundary)),
        "wound_vertices_used": int(len(selected_depth)),
        "projected_opening_area_mm2": _polygon_area(rim_uv),
        "wound_surface_area_mm2": float(mesh.area_faces[selected_faces].sum()),
        "perimeter_mm": perimeter,
        "maximum_length_or_width_mm": _feret_diameter(rim_uv),
        "cavity_volume_mm3": float((face_depth[selected_faces] * projected_face_area[selected_faces]).sum()),
        "maximum_depth_mm": max_depth,
        "mean_depth_mm": float(selected_depth.mean()),
        "median_depth_mm": float(np.median(selected_depth)),
        "depth_p95_mm": float(np.percentile(selected_depth, 95)),
        "depth_standard_deviation_mm": float(selected_depth.std()),
        "mean_wall_slope_degrees": float(slopes.mean()) if len(slopes) else None,
        "maximum_wall_slope_degrees": float(slopes.max()) if len(slopes) else None,
        "reference_surface": "thin-plate interpolation through the traced rim",
        "cavity_direction": "selected automatically from the larger side of the reference surface",
        "width_profile": depth_profile,
    }
