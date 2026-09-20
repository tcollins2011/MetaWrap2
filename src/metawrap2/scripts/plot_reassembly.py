#!/usr/bin/env python
# USAGE:
# ./script file1.stats file2.stats file3.stats

import os
import sys
from typing import List, Tuple

import matplotlib.pyplot as plt

plt.switch_backend("agg")


max_contamination = int(sys.argv[3])
min_completion = int(sys.argv[2])


####################################################################################################################################
############################################         MAKE THE N50 PLOT                 ############################################
####################################################################################################################################
def _bin_set_label(path):
    """Bin-set name from a .stats.tsv path: basename with the suffix removed.

    Both this and the colour assignment below must agree on the label, or a bin set gets no
    colour. Handles the ".stats.tsv" suffix as one unit rather than stripping a single
    extension (which would leave "binsA.stats").
    """
    base = os.path.basename(path)
    for suffix in (".stats.tsv", ".stats"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return os.path.splitext(base)[0]


print("Loading completion info....")
data = {}
max_n50 = 0
# loop over all bin .stats files
for file_name in sys.argv[4:]:
    bin_set = _bin_set_label(file_name)
    data[bin_set] = []
    with open(file_name) as fh:
        lines = fh.readlines()
    for line in lines:
        # skip header
        if "compl" in line:
            continue

        # skip bins that are too contaminated or very incomplete
        if float(line.split("\t")[2]) > max_contamination:
            continue
        if float(line.split("\t")[1]) < min_completion:
            continue
        if float(line.split("\t")[1]) < 1:
            continue

        # save the completion value of each bin into a list
        data[bin_set].append(int(line.split("\t")[5]))
        max_n50 = max(max_n50, int(line.split("\t")[5]))
for values in data.values():
    values.sort(reverse=True)

print("Plotting completion data...")
# MAKING THE PLOT PRETTY!!!!
# set some color schemes
tableau20 = [
    (214, 39, 40),
    (31, 119, 180),
    (255, 127, 14),
    (255, 0, 0),
    (0, 0, 255),
    (197, 176, 213),
    (140, 86, 75),
    (196, 156, 148),
    (227, 119, 194),
    (247, 182, 210),
    (127, 127, 127),
    (199, 199, 199),
    (188, 189, 34),
    (219, 219, 141),
    (23, 190, 207),
    (158, 218, 229),
]

# matplotlib wants 0-1 floats, but the table above is written in 0-255 ints because that is
# how the palette is published; build a new list rather than rewriting it in place so the
# element type stays consistent.
palette: List[Tuple[float, float, float]] = [
    (r / 255.0, g / 255.0, b / 255.0) for r, g, b in tableau20
]
plot_colors = {}
# sys.argv[1] and [2] are the completeness/contamination thresholds, not stats files:
# enumerating argv[1:] shifted every colour by two and, with more than ~15 bin sets, ran off
# the end of the palette with an IndexError. Cycling keeps any number of sets plottable.
for i, label in enumerate(sys.argv[4:]):
    plot_colors[_bin_set_label(label)] = palette[i % len(palette)]


# set figure size
plt.figure(figsize=(18, 8))
plt.style.use("ggplot")

# Remove the plot frame lines. They are unnecessary chartjunk.
ax = plt.subplot(131)
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_linewidth(0.5)
ax.spines["bottom"].set_color("black")
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.set_facecolor("white")

# Ensure that the axis ticks only show up on the bottom and left of the plot.
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

# Limit the range of the plot to only where the data is.
plt.ylim(0, max_n50)
max_x = 0
for values in data.values():
    max_x = max(max_x, len(values))
plt.xlim(0, max_x)

# Make sure your axis ticks are large enough to be easily read.
plt.yticks(
    range(0, max_n50, 50000), [str(x // 1000) + "kb" for x in range(0, max_n50, 50000)], fontsize=14
)
plt.xticks(fontsize=14)

# Provide tick lines across the plot to help your viewers trace along
for y in range(0, max_n50, 50000):
    plt.plot(range(max_x), [y] * len(range(max_x)), "--", lw=0.5, color="black", alpha=0.3)

# Remove the tick marks; they are unnecessary with the tick lines we just plotted.
plt.tick_params(
    axis="both",
    which="both",
    bottom=False,
    top=False,
    labelbottom=True,
    left=False,
    right=False,
    labelleft=True,
)


# PLOTTING THE DATA

# prepare labeles
labels = list(data)


# start plotting data
for rank, bin_set in enumerate(labels):
    # chose a color!
    c = plot_colors[bin_set]

    # plot the data
    plt.plot(data[bin_set], lw=2.5, color=c)

    # add bin set label to plot
    x_pos = len(data[bin_set]) // 4
    y_pos = data[bin_set][x_pos - 1]
    plt.text(x_pos, y_pos, bin_set, fontsize=18, color=c)

# add plot and axis titles and adjust edges
plt.title("Bin N50 ranking", fontsize=26)
plt.xlabel("Descending N50 rank", fontsize=20)
plt.ylabel("Bin N50", fontsize=20)


####################################################################################################################################
############################################         MAKE THE COMPLETION PLOT           ############################################
####################################################################################################################################
print("Loading completion info....")


data = {}
max_x = 0
# loop over all bin .stats files
for file_name in sys.argv[4:]:
    bin_set = _bin_set_label(file_name)
    data[bin_set] = []
    with open(file_name) as fh:
        lines = fh.readlines()
    for line in lines:
        # skip header
        if "compl" in line:
            continue

        # skip bins that are too contaminated or very incomplete
        if float(line.split("\t")[2]) > max_contamination:
            continue
        if float(line.split("\t")[1]) < min_completion:
            continue
        if float(line.split("\t")[1]) < 1:
            continue

        # save the completion value of each bin into a list
        data[bin_set].append(float(line.split("\t")[1]))
    max_x = max(max_x, len(data[bin_set]))

# sort the completion data sets
for bin_set in data:
    data[bin_set].sort(reverse=True)

print("Plotting completion data...")
# set figure size
plt.style.use("ggplot")

# Remove the plot frame lines. They are unnecessary chartjunk.
ax = plt.subplot(132)
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_linewidth(0.5)
ax.spines["bottom"].set_color("black")
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.set_facecolor("white")

# Ensure that the axis ticks only show up on the bottom and left of the plot.
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

# Limit the range of the plot to only where the data is.
plt.ylim(min_completion, 102)
max_x = 0
for k in data:
    max_x = max(max_x, len(data[k]))
plt.xlim(0, max_x)

# Make sure your axis ticks are large enough to be easily read.
plt.yticks(
    range(min_completion, 105, 10),
    [str(x) + "%" for x in range(min_completion, 105, 10)],
    fontsize=14,
)
plt.xticks(fontsize=14)

# Provide tick lines across the plot to help your viewers trace along
for y in range(min_completion, 105, 10):
    plt.plot(range(max_x), [y] * len(range(max_x)), "--", lw=0.5, color="black", alpha=0.3)

# Remove the tick marks; they are unnecessary with the tick lines we just plotted.
plt.tick_params(
    axis="both",
    which="both",
    bottom=False,
    top=False,
    labelbottom=True,
    left=False,
    right=False,
    labelleft=True,
)


# PLOTTING THE DATA

# prepare labeles
labels = []
for k in data:
    labels.append(k)

# make ranking system for lable distribution
ranks = {}
for bin_set in labels:
    p = len(data[bin_set])
    ranks[bin_set] = p
rank_order = {}
n = 0
for n, (key, _value) in enumerate(
    sorted(ranks.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)
):
    rank_order[key] = n


# start plotting data

for rank, bin_set in enumerate(labels):
    # chose a color!
    c = plot_colors[bin_set]

    # plot the data
    plt.plot(data[bin_set], lw=2.5, color=c)

    # add bin set label to plot
    x_pos = len(data[bin_set]) // 2
    if "reasse" in bin_set:
        x_pos = len(data[bin_set]) // 3
    y_pos = data[bin_set][x_pos]
    plt.text(x_pos, y_pos, bin_set, fontsize=18, color=c)

# add plot and axis titles and adjust edges
plt.title("Bin completion ranking", fontsize=26)
plt.xlabel("Descending completion rank", fontsize=20)
plt.ylabel("Estimated bin completion", fontsize=20)


####################################################################################################################################
############################################         MAKE THE CONTAMINATION PLOT        ############################################
####################################################################################################################################
print("Loading contamination info...")

data = {}
# loop over all bin .stats files
for file_name in sys.argv[4:]:
    bin_set = _bin_set_label(file_name)
    data[bin_set] = []
    with open(file_name) as fh:
        lines = fh.readlines()
    for line in lines:
        # skip header
        if "compl" in line:
            continue

        # skip bins that are too incomplete or way too contaminated
        if float(line.split("\t")[1]) < min_completion:
            continue
        if float(line.split("\t")[2]) > max_contamination:
            continue

        # save the contamination value of each bin into a list
        data[bin_set].append(float(line.split("\t")[2]))

# sort the contamination data sets
for bin_set in data:
    data[bin_set].sort(reverse=False)

print("Plotting the contamination data...")
# MAKING THE PLOT PRETTY!!!!
# Remove the plot frame lines. They are unnecessary chartjunk.
ax = plt.subplot(133)
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_linewidth(0.5)
ax.spines["bottom"].set_color("black")
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
ax.set_facecolor("white")

# Ensure that the axis ticks only show up on the bottom and left of the plot.
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

# Limit the range of the plot to only where the data is.
# plt.gca().invert_yaxis()
plt.ylim(0, max_contamination + 0.5)
# ax.set_yscale('log')
max_x = 0
for k in data:
    max_x = max(max_x, len(data[k]))
plt.xlim(0, max_x)

# Make sure your axis ticks are large enough to be easily read.
plt.yticks(
    range(0, max_contamination + 1, 2),
    [str(x) + "%" for x in range(0, max_contamination + 1, 2)],
    fontsize=14,
)
plt.xticks(fontsize=14)

# Provide tick lines across the plot to help your viewers trace along
for y in range(0, max_contamination + 1, 2):
    plt.plot(range(max_x), [y] * len(range(max_x)), "--", lw=0.5, color="black", alpha=0.3)

# Remove the tick marks; they are unnecessary with the tick lines we just plotted.
plt.tick_params(
    axis="both",
    which="both",
    bottom=False,
    top=False,
    labelbottom=True,
    left=False,
    right=False,
    labelleft=False,
)


# PLOTTING THE DATA
# prepare labeles
labels = []
for k in data:
    labels.append(k)


# start plotting data
for rank, bin_set in enumerate(labels):
    # chose a color!
    c = plot_colors[bin_set]

    # plot the data
    plt.plot(data[bin_set], lw=2.5, color=c)

    # add plot label
    # x_pos = len(data[bin_set])-1-20*(len(rank_order)-rank_order[bin_set]-1)
    x_pos = len(data[bin_set]) // 3
    y_pos = data[bin_set][x_pos]
    plt.text(x_pos + 1, y_pos, bin_set, fontsize=18, color=c)

# add plot and axis titles and adjust the edges
plt.title("Bin contamination ranking", fontsize=26)
plt.xlabel("Acending contamination rank", fontsize=20)
plt.ylabel("Estimated bin contamination", fontsize=20)
plt.gcf().subplots_adjust(right=0.9)


# save figure
print("Saving figures reassembly_results.eps and reassembly_results.png to folder " + sys.argv[1])
plt.tight_layout(w_pad=5)
# plt.subplots_adjust(top=0.92, right=0.90, left=0.08)
plt.savefig(sys.argv[1] + "/" + "reassembly_results.eps", format="eps", dpi=600)
plt.savefig(sys.argv[1] + "/" + "reassembly_results.png", format="png", dpi=600)
# plt.show()
