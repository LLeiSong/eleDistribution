from pylandstats import *


# Initialize landscape
fname = 'results/tanzania/prediction/landcover.tif'
landscape = Landscape(fname)

# Patch-level metrics
patch_areas = landscape.area()
patch_perimeter = landscape.perimeter()
patch_par = landscape.perimeter_area_ratio()
patch_shape_index = landscape.shape_index()
patch_fractal = landscape.fractal_dimension()
patch_nn = landscape.euclidean_nearest_neighbor()

patch_metrics = landscape.compute_patch_metrics_df()

# Class-level metrics
class_area = landscape.total_area()
class_proportion = landscape.proportion_of_landscape()
class_num_patches = landscape.number_of_patches()
class_patch_density = landscape.patch_density()
class_largest_patch = landscape.largest_patch_index()
class_total_edge = landscape.total_edge()
class_edge_density = landscape.edge_density()
class_shape_index = landscape.landscape_shape_index()
class_mesh_size = landscape.effective_mesh_size()

class_metrics = landscape.compute_class_metrics_df()

# Landscape-level metrics
landscape_entropy = landscape.entropy()
landscape_diversity_index = landscape.shannon_diversity_index()
landscape_joint_entropy = landscape.joint_entropy()
landscape_conditional_entropy = landscape.conditional_entropy()
landscape_mutual_info = landscape.mutual_information()
landscape_relative_mutual_info = landscape.relative_mutual_information()
landscape_contagion = landscape.contagion()

landscape_metrics = landscape.compute_landscape_metrics_df()
