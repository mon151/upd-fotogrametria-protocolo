# upd-fotogrametria-protocolo
Pipeline de análisis geométrico para el protocolo de fotogrametría de úlceras de pie diabético | Proyecto Integrador II
## Módulos
- `config.py` — parámetros del sistema (archivo del modelo, tamaño de vóxel, batch)
- `main.py` — orquesta el cálculo de volumen por voxelización
- `voxelizer.py` — clase `VolumeCalculator`: discretiza la malla en vóxeles y cuenta los centros contenidos dentro del sólido
- `wound_boundary.py` — plantilla CSV para el borde de la herida trazado manualmente
- `wound_metrics.py` — cálculo de área, profundidad, volumen y perfil de ancho sobre la región segmentada
- `main_wound_metrics.py` — script ejecutable que corre el análisis y guarda un reporte JSON
- `concavity_main.py` — cálculo del volumen de concavidad/socavamiento
- ## Requisitos
- python -m pip install numpy trimesh rtree pyvista vtk scipy meshio matplotlib pandas networkx tqdm open3d scikit-image opencv-python
