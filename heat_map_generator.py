"""
Heat Map Generator by Eric Meunier

This module generates a heat map of the proximity of any two rock types.
This is accomplished by creating a grid of coordinates over the zone of interest, and calculating the distance to the nearest two selected rock units. 
The mean distance between the two will be used as the distance value for this coordinate. This mesh will then be used to create a contour map.
Distance to a geological unit is calculated by calculating the distance from each grid mesh point to each coordinate that defines the exterior of the polygon and using the minima.

Since there are many calculations to be made, vectorized numpy functions are used wherever possible to save on computation time.

Requires geopandas and matplotlib.

Installation of dependencies:
    pip install pyogrio
    pip install geopandas
    pip install matplotlib

Usage:
    Create an instance of HeatMapGenerator, and call the 'generate_heatmap' method by passing in the shapefile filepath (as a string or pathlib.Path), 
    the details as a dictionary (column name in the shapefile as the key, what string to match as the value) of two rock types, and an optional falloff distance (as an integer). 
    If no falloff distance is specified, 10,000 will be used.
    
    The distance value beneath the mouse cursor is shown in the bottom right of the window. Additional, clicking on the map will also print the z value to the console.

    Example: 
        shapefile = r"BedrockP.shp"

        rock_1_details = {"rock_type": "grano"}
        rock_2_details = {"rock_type": "serp", "rock_class": "ultramafic rocks"}
        
        map = HeatMapGenerator()
        map.generate_heatmap(shapefile, rock_1_details, rock_2_details, falloff_distance=8000)
        map.show()
        map.save("Colbat heat map.png")
"""

import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from mpl_toolkits.axes_grid1 import make_axes_locatable


def find_nearest(target_value: float or int, arr: np.array) -> float or int:
    """
    Return the item in an array closest in value to the target value.
    """
    arr = np.asarray(arr)
    idx = (np.abs(arr - target_value)).argmin()
    return arr[idx]


def calc_distance(target_x: float, target_y: float, xx: np.array, yy: np.array) -> np.array:
    """
    Vectorized calculation of the euclidean distance between a target coordinate and each coordinate in the xx, yy grid mesh
    """
    return np.sqrt((xx - target_x) ** 2 + (yy - target_y) ** 2)


def get_all_coordinates(df: gpd.GeoDataFrame) -> tuple:
    """
    Retrieve all x, y coordinates of the exterior of every geometry in a geodataframe.
    """
    xs, ys = [], []
    for geom in df.geometry:
        x, y = geom.exterior.coords.xy
        xs.extend(x)
        ys.extend(y)
    return xs, ys


def get_min_distance(xs: list, ys: list, xx: np.array, yy: np.array) -> np.array:
    """
    Returns the set of distances where the smallest distance to a grid point was found.
    """
    zs = []
    for x, y in zip(xs, ys):
        z = calc_distance(x, y, xx, yy)
        zs.append(z)
    return np.array(zs).min(axis=0)


def create_grid_mesh(df: gpd.GeoDataFrame, falloff_distance: int=0) -> tuple:
    """
    Creates all combination of X and Y coordinates to form a 2D grid covering the geometry of a geodataframe (including any falloff distance).
    """
    # Apply a buffer to the zone of interest to ensure our measurements go at least to the falloff distance.
    xmin, ymin, xmax, ymax = df.buffer(falloff_distance).total_bounds
    print(f"Grid bounds: {xmin:.0f}, {ymin:.0f}, {ymin:.0f}, {ymax:.0f}")
    grid_spacing = 2500  # meters

    xi = np.arange(xmin, xmax + grid_spacing, grid_spacing)  # Add grid_spacing to the maximum value because np.arange excludes the upper limit
    yi = np.arange(ymin, ymax + grid_spacing, grid_spacing)
    xx, yy = np.meshgrid(xi, yi)
    return xi, yi, xx, yy
    

class HeatMapGenerator:
    def __init__(self):
        """
        Class which generates a heat map based on the proximity of two rock types from an input shapefile.
        """
        self.falloff_distance = None

        self.fig, self.ax = plt.subplots(figsize=(10, 10))
        self.ax.set_xlabel("Easting (m)")
        self.ax.set_ylabel("Northing (m)")
        self.ax.yaxis.get_major_formatter().set_scientific(False)  # Disable scientific notation for the Northings
        self.ax.format_coord = self._format_coord    
        divider = make_axes_locatable(self.ax)
        self.cax = divider.append_axes("right", size="5%", pad=0.2)
        self.cmap = plt.get_cmap('jet').reversed()
        self.cbar = None
        plt.connect('button_press_event', self._on_click)
            
    def _calculate_distances(self, df1: gpd.GeoDataFrame, df2: gpd.GeoDataFrame) -> np.array:
        df1_xs, df1_ys = get_all_coordinates(df1)    
        df2_xs, df2_ys = get_all_coordinates(df2) 

        df1_zz = get_min_distance(df1_xs, df1_ys, self.xx, self.yy)
        df2_zz = get_min_distance(df2_xs, df2_ys, self.xx, self.yy)

        return (df1_zz + df2_zz) / 2  # Mean distance

    def _format_coord(self, x: float, y: float):
        """
        Overload the function for displaying the mouse coordinates when mousing-over in the matplotlib plot to include the z value.
        """
        z = self.get_z(x, y)
        return f"x= {x:.2f}, y= {y:.2f}, z= {z:.2f}"

    def _on_click(self, event):
        """
        Print the z value of the coordinate nearest the click event.
        """
        if event.inaxes:
            z = self.get_z(event.xdata, event.ydata)
            print(f"{z:.2f}")

    def _plot_heatmap(self):
        """
        Create the figure and plot the shapefile units along with the contour map of the distances. Nothing is plotted beyond the falloff distance.
        """
        norm = mpl.colors.Normalize(vmin=0, vmax=self.falloff_distance)
        self.cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=self.cmap), ax=self.ax, cax=self.cax)
        self.cbar.ax.invert_yaxis()  # Invert the colorbar since the smaller value equals a higher 

        self.ax.contourf(self.xx, self.yy, self.zz, levels=np.linspace(0, self.falloff_distance, 100), alpha=0.2, cmap=self.cmap, zorder=5)
        self.shp_df.plot(ax=self.ax, column="era", categorical=True, legend=False, cmap='copper')

    def get_z(self, x: float, y: float) -> float:
        """
        Return the z value of the coordinate nearest x and y.
        """
        coord_index = np.argwhere((self.xx==find_nearest(x, self.xi)) & (self.yy==find_nearest(y, self.yi)))
        z = self.zz[coord_index[0][0], coord_index[0][1]]
        return z

    def set_title(self, title: str):
        self.ax.set_title(title)

    def set_cbar_label(self, label: str):
        self.cbar.set_label(label)

    def show(self):
        plt.show()

    def save(self, save_name: str):
        self.fig.savefig(save_name)

    def generate_heatmap(self, shp_file: str, rock_1_details: list, rock_2_details: list, falloff_distance: int=10000):
        print(f"Generating heatmap for {Path(shp_file).name}. Falloff distance: {falloff_distance}.")
        self.falloff_distance = falloff_distance
        self.shp_df = gpd.read_file(shp_file)

        rock_1_filter_str = ' | '.join([f"{column}.str.contains('{value}')" for (column, value) in rock_1_details.items()])
        rock_2_filter_str = ' | '.join([f"{column}.str.contains('{value}')" for (column, value) in rock_2_details.items()])
        rock_1_df = self.shp_df.query(rock_1_filter_str)
        rock_2_df = self.shp_df.query(rock_2_filter_str)
        assert not rock_1_df.empty, f"No entries were found for the query '{rock_1_filter_str}'"
        assert not rock_2_df.empty, f"No entries were found for the query '{rock_2_filter_str}'"
        
        # Only calculate distances that are within the zone of interest (including the falloff distance) to reduce computation time.
        zone_of_interest = self.shp_df.query(f"{rock_1_filter_str} or {rock_2_filter_str}")  
        self.xi, self.yi, self.xx, self.yy = create_grid_mesh(zone_of_interest, self.falloff_distance)
        self.zz = self._calculate_distances(rock_1_df, rock_2_df)
        self._plot_heatmap()


if __name__ == "__main__":
    shapefile = r"Sample_files/BedrockP.shp"

    rock_1_details = {"rock_type": "grano"}
    rock_2_details = {"rock_type": "serp", "rock_class": "ultramafic rocks"}
    
    map = HeatMapGenerator()
    map.generate_heatmap(shapefile, rock_1_details, rock_2_details, falloff_distance=8000)
    map.set_title("Probability of finding Cobalt")
    map.set_cbar_label("Mean distance from Granodiorite and Serpentinite (m)\nSmaller value equates to higher probability")
    map.show()
    map.save(f"Cobalt heat map.png")
