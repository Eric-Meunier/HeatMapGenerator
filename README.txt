Heat Map Generator by Eric Meunier

There are two python script files here, one is my working copy which shows step-by-step how I 
solved the problem with some commentary (Cobalt heat map (working copy).ipynb, along with an HTML version so it can be viewed more easily). 
The other (heat_map_generator.py) is the more final deliverable which includes the heatmap generator tool. 
A static heat map image is also included.

The generator tool creates a heat map of the proximity of any two rock types.
This is accomplished by creating a grid of coordinates over the zone of interest, and calculating the 
distance to the nearest two selected rock units. 
The mean distance between the two will be used as the distance value for this coordinate. 
This mesh will then be used to create a contour map.
Distance to a geological unit is calculated by calculating the distance from each grid mesh point 
to each coordinate that defines the exterior of the polygon and using the minima.

Since there are many calculations to be made, vectorized numpy functions are used wherever possible 
to save on computation time.

Requires geopandas and matplotlib.

Installation of dependencies:
    pip install pyogrio
    pip install geopandas
    pip install matplotlib

Usage:
    Create an instance of HeatMapGenerator, and call the 'generate_heatmap' method by passing in the 
	shapefile filepath (as a string or pathlib.Path), the details as a dictionary (column name in 
	the shapefile as the key, what string to match as the value) of two rock types, and an optional 
	falloff distance (as an integer). If no falloff distance is specified, 10,000 will be used.

    The distance value beneath the mouse cursor is shown in the bottom right of the window. Additional, clicking on the map will also print the z value to the console.

    Example: 
        shapefile = r"BedrockP.shp"

        rock_1_details = {"rock_type": "grano"}
        rock_2_details = {"rock_type": "serp", "rock_class": "ultramafic rocks"}
        
        map = HeatMapGenerator()
        map.generate_heatmap(shapefile, rock_1_details, rock_2_details, falloff_distance=8000)
        map.show()