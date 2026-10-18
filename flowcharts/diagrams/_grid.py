"""One grid for every diagram, so boxes are the same size and the same distance apart."""

W = 290  # step box width
COL = 340  # column pitch: 50 px between boxes
ROW = 135  # row pitch


def x(col: int) -> float:
    return 170 + COL * col


def y(row: float) -> float:
    return 120 + ROW * row


#: Panel colours (fill, stroke), from the module palette lightened.
RED = ("#FDF1F0", "#B85450")
BLUE = ("#F0F6FC", "#6C8EBF")
GREEN = ("#F2F8EF", "#82B366")
YELLOW = ("#FFF9E6", "#D6B656")
PURPLE = ("#F7F1FA", "#9673A6")
