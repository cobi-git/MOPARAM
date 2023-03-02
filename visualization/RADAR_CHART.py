import pathlib
from os.path import join

import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Circle, RegularPolygon
from matplotlib.path import Path
from matplotlib.projections.polar import PolarAxes
from matplotlib.projections import register_projection
from matplotlib.spines import Spine
from matplotlib.transforms import Affine2D

from library.TEST_UTILS import TEST_OMICS, CLUSTERING_SCORE_DIR, OUT_BASE_DIR
from library.modb_utils import getdirlist
from visualization.print_test_result import data_load


def radar_factory(num_vars, frame='circle'):
    # calculate evenly-spaced axis angles
    theta = np.linspace(0, 2*np.pi, num_vars, endpoint=False)

    class RadarTransform(PolarAxes.PolarTransform):

        def transform_path_non_affine(self, path):
            # Paths with non-unit interpolation steps correspond to gridlines,
            # in which case we force interpolation (to defeat PolarTransform's
            # autoconversion to circular arcs).
            if path._interpolation_steps > 1:
                path = path.interpolated(num_vars)
            return Path(self.transform(path.vertices), path.codes)

    class RadarAxes(PolarAxes):

        name = 'radar'
        # use 1 line segment to connect specified points
        RESOLUTION = 1
        PolarTransform = RadarTransform

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            # rotate plot such that the first axis is at the top
            self.set_theta_zero_location('N')

        def fill(self, *args, closed=True, **kwargs):
            """Override fill so that line is closed by default"""
            return super().fill(closed=closed, *args, **kwargs)

        def plot(self, *args, **kwargs):
            """Override plot so that line is closed by default"""
            lines = super().plot(*args, **kwargs)
            for line in lines:
                self._close_line(line)

        def _close_line(self, line):
            x, y = line.get_data()
            # FIXME: markers at x[0], y[0] get doubled-up
            if x[0] != x[-1]:
                x = np.append(x, x[0])
                y = np.append(y, y[0])
                line.set_data(x, y)

        def set_varlabels(self, labels):
            self.set_thetagrids(np.degrees(theta), labels)

        def _gen_axes_patch(self):
            # The Axes patch must be centered at (0.5, 0.5) and of radius 0.5
            # in axes coordinates.
            if frame == 'circle':
                return Circle((0.5, 0.5), 0.5)
            elif frame == 'polygon':
                return RegularPolygon((0.5, 0.5), num_vars,
                                      radius=.5, edgecolor="k")
            else:
                raise ValueError("Unknown value for 'frame': %s" % frame)

        def _gen_axes_spines(self):
            if frame == 'circle':
                return super()._gen_axes_spines()
            elif frame == 'polygon':
                # spine_type must be 'left'/'right'/'top'/'bottom'/'circle'.
                spine = Spine(axes=self,
                              spine_type='circle',
                              path=Path.unit_regular_polygon(num_vars))
                # unit_regular_polygon gives a polygon of radius 1 centered at
                # (0, 0) but we want a polygon of radius 0.5 centered at (0.5,
                # 0.5) in axes coordinates.
                spine.set_transform(Affine2D().scale(.5).translate(.5, .5)
                                    + self.transAxes)
                return {'polar': spine}
            else:
                raise ValueError("Unknown value for 'frame': %s" % frame)

    register_projection(RadarAxes)
    return theta


def load_data(test):
    if test == 'omics_combination_test':
        total = []
        tool = None
        cancer_list = getdirlist(join(score_path, test))
        df_all = pd.DataFrame()
        for i, cancer in enumerate(cancer_list):
            cluster_score_path = join(score_path, test, cancer, "cluster_score.norm.csv")
            df = pd.read_csv(cluster_score_path, index_col=0).drop(columns=['Best', 'Mean'])
            df = df.transpose()
            #print(df)
            if i == 0:
                total.append(df.columns)
                tool = df.index
                df_all = df
            else:
                df_all = df_all + df
            total.append((cancer, df.values.tolist()))
        df_all = df_all/len(cancer_list)
        total.insert(1, ('Total mean', df_all.mean().tolist()))
        return total, len(total[0]), tool
    else:
        tool_dict = dict()
        cancer_list = getdirlist(join(score_path, test))
        total = []
        df_all = pd.DataFrame()
        for i, cancer in enumerate(cancer_list):
            cluster_score_path = join(score_path, test, cancer, "cluster_score.norm.csv")
            df = pd.read_csv(cluster_score_path, index_col=0).drop(columns=['Best', 'Mean'])
            if i == 0:
                total.append(df.index.tolist())

            for tool in df.columns:
                if not tool in tool_dict:
                    tool_dict[tool] = pd.DataFrame()
                tool_dict[tool][cancer] = df[tool]

        for tool, df in tool_dict.items():
            df = df.transpose()
            total.append((tool, df.values.tolist()))
        df_all = sum(tool_dict.values())/len(tool_dict)
        df_all = df_all.transpose()
        total.insert(1, ('Total mean', df_all.mean().tolist()))
        return total, len(total[0]), cancer_list

if __name__ == '__main__':
    score_path = CLUSTERING_SCORE_DIR
    test_list = getdirlist(score_path)
    plotdir = join(OUT_BASE_DIR, 'Radar_chart')
    pathlib.Path(plotdir).mkdir(parents=True, exist_ok=True)
    for test in test_list:
        if test.find('omics_combination') == -1:# only for omics combination
            continue

        print(test)
        data, N, tool = load_data(test)
        #print(data)
        theta = radar_factory(N, frame='polygon')
        spoke_labels = data.pop(0)

        #fig, axs = plt.subplots(figsize=(9, 18), nrows=6, ncols=2,
        #                        subplot_kw=dict(projection='radar'))

        fig = plt.figure(figsize=(9,18))
        gs = fig.add_gridspec(6,4)
        fig.subplots_adjust(wspace=0.25, hspace=0.20, top=0.85, bottom=0.05)
        colors = ['C' + str(a) for a in range(0,10,1)]
        # Plot the four cases from the example data on separate axes
        for i, (title, case_data) in enumerate(data):
            col = 2*(i%2)
            if i == 0:
                ax = plt.subplot(gs[0, col+1:col+3], projection='radar')#total mean
                max_point = round(max(case_data), 1) + 0.1
            else:
                ax = plt.subplot(gs[int((i + 1)/2), col:col + 2], projection='radar')#else
                max_point = round(max([max(a) for a in case_data]), 1) + 0.1
            ax.set_rgrids([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

            ax.set_ylim(0.0, max_point)
            ax.set_title(title, weight='bold', size='medium', position=(0.5, 1.1),
                         horizontalalignment='center', verticalalignment='center')

            if i == 0:
                ax.plot(theta, case_data, color='r')
                ax.fill(theta, case_data, facecolor='r', alpha=0.25, label='_nolegend_')
            else:
                for d, color in zip(case_data, colors):
                    ax.plot(theta, d, color=color)
                    ax.fill(theta, d, facecolor=color, alpha=0.25, label='_nolegend_')

            ax.set_varlabels(spoke_labels)
        # if test == 'omics_combination_test_gender':
        #     fig.delaxes(fig.axes[-1])

        # add legend relative to top-left plot
        labels = tool
        legend = fig.legend(labels, loc="lower center", labelspacing=0.1, fontsize='small', mode='expand', ncol=5)
        #fig.text(0.5, 0.965, test, horizontalalignment='center', color='black', weight='bold', size='large')
        fig.tight_layout()
        plt.subplots_adjust(bottom=0.03)
        plt.savefig(join(plotdir, test + '.radar.png'))