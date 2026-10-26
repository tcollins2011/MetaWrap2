#!/usr/bin/env python
# USAGE:
# ./script MIN_COMPLETION MAX_CONTAMINATION file1.stats file2.stats file3.stats
# Writes binning_results.png in the current directory.

import os
import sys
from typing import List, Tuple

import matplotlib.pyplot as plt

plt.switch_backend("agg")

max_contamination = int(sys.argv[2])
min_completion = int(sys.argv[1])

####################################################################################################################################
############################################         MAKE THE COMPLETION PLOT           ############################################
####################################################################################################################################
print("Loading completion info....")


def _bin_set_label(path):
    """Bin-set name from a .stats.tsv path: basename with the suffix removed.

    Also used for the colour assignment below, so both must give the same label. ".stats.tsv" is
    removed as one suffix, so "binsA.stats.tsv" gives "binsA".
    """
    base = os.path.basename(path)
    for suffix in (".stats.tsv", ".stats"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return os.path.splitext(base)[0]


data = {}
max_x = 0
# loop over all bin .stats files
for file_name in sys.argv[3:]:
    print(file_name)
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

        # save the completion value of each bin into a list
        data[bin_set].append(float(line.split("\t")[1]))
    max_x = max(max_x, len(data[bin_set]))

# sort the completion data sets
for values in data.values():
    values.sort(reverse=True)

print("Plotting completion data...")
# Plot styling.
# set some color schemes
tableau20 = [
    (214, 39, 40),
    (31, 119, 180),
    (255, 127, 14),
    (44, 160, 44),
    (255, 152, 150),
    (148, 103, 189),
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

# Convert the published 0-255 palette to the 0-1 floats matplotlib expects.
palette: List[Tuple[float, float, float]] = [
    (r / 255.0, g / 255.0, b / 255.0) for r, g, b in tableau20
]
plot_colors = {}
# Stats files start at argv[3] (argv[1] and [2] are thresholds). The palette is cycled so any
# number of bin sets can be plotted.
for i, label in enumerate(sys.argv[3:]):
    plot_colors[_bin_set_label(label)] = palette[i % len(palette)]


# set figure size
plt.figure(figsize=(16, 8))
plt.style.use("ggplot")

# Remove the top, right and left frame lines.
ax = plt.subplot(121)
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_linewidth(0.5)
ax.spines["bottom"].set_color("black")
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
# ax.set_facecolor('white')
ax.set_facecolor("white")

# Ensure that the axis ticks only show up on the bottom and left of the plot.
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

# Limit the range of the plot to only where the data is.
plt.ylim(min_completion, 105)
max_x = 0
for values in data.values():
    max_x = max(max_x, len(values))
plt.xlim(0, max_x)

# Tick label size.
plt.yticks(
    range(min_completion, 105, 10),
    [str(x) + "%" for x in range(min_completion, 105, 10)],
    fontsize=14,
)
plt.xticks(fontsize=14)

# Grid lines.
for y in range(min_completion, 105, 10):
    plt.axhline(y=y, linestyle="--", lw=0.5, color="black", alpha=0.3)
for x in range(0, 1000, 20):
    plt.axvline(x=x, linestyle="--", lw=0.5, color="black", alpha=0.3)

# Hide tick marks; the grid lines replace them.
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

# prepare labels
labels = list(data)

# plot the data and labels
N = len(labels)
y_increment = (100 - min_completion) // N // 2
y_pos = 100 - y_increment

for rank, bin_set in enumerate(labels):
    # bin set colour
    c = plot_colors[bin_set]

    # plot the data
    plt.plot(data[bin_set], lw=2.5, color=c)

    # add bin set label to plot; an empty bin set has no curve, so its label goes at x = 0
    x_pos = 0
    for x_pos, y in enumerate(data[bin_set]):
        if y < y_pos:
            break
    plt.text(x_pos, y_pos, bin_set, fontsize=18, color=c)
    y_pos -= y_increment

# add plot and axis titles and adjust edges
plt.title("Bin completion ranking", fontsize=26)
plt.xlabel("Descending completion rank", fontsize=16)
plt.ylabel("Estimated bin completion", fontsize=16)


####################################################################################################################################
############################################         MAKE THE CONTAMINATION PLOT        ############################################
####################################################################################################################################
print("Loading contamination info...")

data = {}
# loop over all bin .stats files
for file_name in sys.argv[3:]:
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
# Plot styling.
# Remove the top, right and left frame lines.
ax = plt.subplot(122)
ax.spines["top"].set_visible(False)
ax.spines["bottom"].set_linewidth(0.5)
ax.spines["bottom"].set_color("black")
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
# ax.set_facecolor('white')
ax.set_facecolor("white")

# Ensure that the axis ticks only show up on the bottom and left of the plot.
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

# Limit the range of the plot to only where the data is.
# plt.gca().invert_yaxis()
plt.ylim(0, max_contamination + 1)
# ax.set_yscale('log')
max_x = 0
for k in data:
    max_x = max(max_x, len(data[k]))
plt.xlim(0, max_x)

# Tick label size.
plt.yticks(
    range(-0, max_contamination + 1, 1),
    [str(x) + "%" for x in range(-0, max_contamination + 1, 1)],
    fontsize=14,
)
plt.xticks(fontsize=14)

# Grid lines.
for y in range(0, max_contamination + 1, 1):
    plt.axhline(y=y, linestyle="--", lw=0.5, color="black", alpha=0.3)
for x in range(0, 1000, 20):
    plt.axvline(x=x, linestyle="--", lw=0.5, color="black", alpha=0.3)


# Hide tick marks; the grid lines replace them.
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
# prepare labels
labels = []
for k in data:
    labels.append(k)

# plot the data and labels
N = len(labels)
y_increment = max_contamination // N // 2
y_pos = y_increment

for rank, bin_set in enumerate(labels):
    # bin set colour
    c = plot_colors[bin_set]

    # plot the data
    plt.plot(data[bin_set], lw=2.5, color=c)

    # add bin set label to plot; an empty bin set has no curve, so its label goes at x = 0
    x_pos = 0
    for x_pos, y in enumerate(data[bin_set]):
        if y > y_pos:
            break
    plt.text(x_pos, y_pos, bin_set, fontsize=18, color=c)
    y_pos += y_increment


# add plot and axis titles and adjust the edges
plt.title("Bin contamination ranking", fontsize=26)
plt.xlabel("Ascending contamination rank", fontsize=16)
plt.ylabel("Estimated bin contamination (log scale)", fontsize=16)
plt.gcf().subplots_adjust(right=0.9)

# save figure
print("Saving figures binning_results.eps and binning_results.png ...")
plt.tight_layout(w_pad=10)
plt.subplots_adjust(top=0.92, right=0.90, left=0.08)
plt.savefig("binning_results.png", format="png", dpi=300)
# plt.show()
