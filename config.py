"""Project settings for the voxel-volume calculation."""

# Put the 3D model in this same folder. OBJ and STL files are supported.
MODEL_FILE = "Model.stl"

# Side length of every voxel, in millimeters.
# 1.0 mm is a good first run. Smaller values improve precision but take longer.
VOXEL_SIZE_MM = 1.0

# Points evaluated at once. This controls memory use and how often the progress
# bar refreshes; it does not change the calculated volume.
POINTS_PER_BATCH = 20_000
