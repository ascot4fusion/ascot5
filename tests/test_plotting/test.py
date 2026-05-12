import pytest
import numpy as np

import a5py.plotting as a5plt


def generate_data(n, mode="linear_pos", seed=0):
    rng = np.random.default_rng(seed)

    log_min, log_max = -3, 3
    if mode == "linear_pos":
        x = rng.random(n)
    elif mode == "linear_signed":
        x = rng.standard_normal(n)
    elif mode == "log_pos":
        x = 10 ** rng.uniform(log_min, log_max, n)
    elif mode == "log_signed":
        mag = 10 ** rng.uniform(log_min, log_max, n)
        sign = rng.choice([-1, 1], size=n)
        x = sign * mag
    else:
        raise ValueError("Unknown mode")
    return x

@pytest.mark.mpl_image_compare
@pytest.mark.parametrize(
    "params",
    [
        {"color": None},
        {"color": "red"},
        {},
        {"scale": "linear_signed"},
        {"scale": "log_pos"},
        {"scale": "log_signed"},
        {"color_interval": np.linspace(-10, 10, 11)},
    ],
    ids=["no-color", "uniform-color", "linear-pos", "linear-signed",
         "log-pos", "log-signed", "custominterval"],
    )
def test_scatter2d(params):
    fig, ax = a5plt.subplots()
    scale = "linear_pos" if "scale" not in params else params["scale"]
    c = generate_data(10, mode=scale) if "color" not in params else params["color"]
    x, y = (
        generate_data(10, mode=scale, seed=0),
        generate_data(10, mode=scale, seed=1),
    )
    color_intervals = 9 if "color_interval" not in params else params["color_interval"]
    scale = "linear" if "linear" in scale else "log"
    a5plt.scatter2d(
        x, y, c, xlabel="x-label", ylabel="y-label", clabel="c-label",
        xscale=scale, yscale=scale, cscale=scale,
        color_intervals=color_intervals, axes=ax,
        )
    return fig

@pytest.mark.mpl_image_compare
@pytest.mark.parametrize(
    "params",
    [
        {"color": None},
        {"color": "red"},
        {},
        {"scale": "linear_signed"},
        {"scale": "log_pos"},
        {"scale": "log_signed"},
        {"color_interval": np.linspace(-10, 10, 11)},
    ],
    ids=["no-color", "uniform-color", "linear-pos", "linear-signed",
         "log-pos", "log-signed", "custominterval"],
    )
def test_scatter3d(params):
    fig, ax = a5plt.subplots(subplot_kw=dict(projection="3d"))
    scale = "linear_pos" if "scale" not in params else params["scale"]
    c = generate_data(10, mode=scale) if "color" not in params else params["color"]
    x, y, z = (
        generate_data(10, mode=scale, seed=0),
        generate_data(10, mode=scale, seed=1),
        generate_data(10, mode=scale, seed=2),
    )
    color_intervals = 9 if "color_interval" not in params else params["color_interval"]
    scale = "linear" if "linear" in scale else "log"
    a5plt.scatter3d(
        x, y, z, c, xlabel="x-label", ylabel="y-label", zlabel="z-label",
        clabel="c-label", xscale=scale, yscale=scale, zscale=scale,
        cscale=scale, color_intervals=color_intervals, axes=ax,
        )
    return fig

@pytest.mark.mpl_image_compare
@pytest.mark.parametrize(
    "params",
    [
        {"color": None},
        {"color": "red"},
        {},
        {"scale": "linear_signed"},
        {"scale": "log_pos"},
        {"scale": "log_signed"},
    ],
    ids=["no-color", "uniform-color", "linear-pos", "linear-signed",
         "log-pos", "log-signed"],
    )
def test_line2d(params):
    scale = "linear_pos" if "scale" not in params else params["scale"]
    c = [generate_data(10, mode=scale)] if "color" not in params else params["color"]
    x, y = (
        generate_data(10, mode=scale, seed=0),
        generate_data(10, mode=scale, seed=1),
    )
    scale = "linear" if "linear" in scale else "log"

    fig, ax = a5plt.subplots()
    a5plt.line2d(
        [x], [y], c, xlabel="x-label", ylabel="y-label",
        clabel="c-label", axes=ax,
        )
    return fig

@pytest.mark.mpl_image_compare
@pytest.mark.parametrize(
    "params",
    [
        {"color": None},
        {"color": "red"},
        {},
        {"scale": "linear_signed"},
        {"scale": "log_pos"},
        {"scale": "log_signed"},
    ],
    ids=["no-color", "uniform-color", "linear-pos", "linear-signed",
         "log-pos", "log-signed"],
    )
def test_line3d(params):
    scale = "linear_pos" if "scale" not in params else params["scale"]
    c = [generate_data(10, mode=scale)] if "color" not in params else params["color"]
    x, y, z = (
        generate_data(10, mode=scale, seed=0),
        generate_data(10, mode=scale, seed=1),
        generate_data(10, mode=scale, seed=2),
    )
    scale = "linear" if "linear" in scale else "log"

    fig, ax = a5plt.subplots(subplot_kw=dict(projection="3d"))
    a5plt.line3d(
        [x], [y], [z], c, xlabel="x-label", ylabel="y-label", zlabel="z-label",
        clabel="c-label", axes=ax,
        )
    return fig

@pytest.mark.mpl_image_compare
@pytest.mark.parametrize(
    "scale",
    ["linear_pos", "linear_signed", "log_pos", "log_signed",],
    ids=["linear-pos", "linear-signed", "log-pos", "log-signed"],
    )
def test_mesh1d(scale):
    x = np.sort(generate_data(20, mode=scale, seed=0))
    y = generate_data(19, mode=scale, seed=1)
    scale = "linear" if "linear" in scale else "log"

    fig, ax = a5plt.subplots()
    a5plt.mesh1d(
        x, y, xscale=scale, yscale=scale, xlabel="x-label", ylabel="y-label",
        axes=ax,
        )
    return fig

@pytest.mark.mpl_image_compare
@pytest.mark.parametrize(
    "scale",
    ["linear_pos", "linear_signed", "log_pos", "log_signed",],
    ids=["linear-pos", "linear-signed", "log-pos", "log-signed"],
    )
def test_mesh2d(scale):
    nx, ny = 20, 10
    x = np.sort(generate_data(nx, mode=scale, seed=0))
    y = np.sort(generate_data(ny, mode=scale, seed=1))
    c1 = np.reshape(
        generate_data(nx*ny, mode=scale, seed=2), (nx, ny),
        )
    nx -= 1
    ny -= 1
    c2 = np.reshape(
        generate_data(nx*ny, mode=scale, seed=2), (nx, ny),
        )
    c1[5, 6] = np.nan
    c2[5, 6] = np.nan
    scale = "linear" if "linear" in scale else "log"
    fig = a5plt.figure()
    ax = fig.add_subplot(1, 2, 1)
    a5plt.mesh2d(
        x, y, c1, xscale=scale, yscale=scale, cscale=scale, xlabel="x-label",
        ylabel="y-label", clabel="c-label", axes=ax,
        )
    ax = fig.add_subplot(1, 2, 2)
    a5plt.mesh2d(
        x, y, c2, xscale=scale, yscale=scale, cscale=scale, xlabel="x-label",
        ylabel="y-label", clabel="c-label", axes=ax,
        )
    return fig

@pytest.mark.mpl_image_compare
@pytest.mark.parametrize(
    "scale",
    ["linear_pos", "linear_signed", "log_pos", "log_signed",],
    ids=["linear-pos", "linear-signed", "log-pos", "log-signed"],
    )
def test_hist1d(scale):
    xbins = np.sort(generate_data(10, mode=scale, seed=0))
    xdata1 = generate_data(1000, mode=scale, seed=1)
    xdata2 = generate_data(500, mode=scale, seed=2)
    wdata1 = generate_data(1000, mode=scale, seed=1)
    wdata2 = generate_data(500, mode=scale, seed=2)
    scale = "linear" if "linear" in scale else "log"

    fig = a5plt.figure()
    ax = fig.add_subplot(1, 2, 1)
    a5plt.hist1d(
        [xdata1, xdata2], xbins=xbins, legend=["1", "2"], xscale=scale,
        xlabel="x-label", ylabel="y-label", axes=ax,
        )
    ax = fig.add_subplot(1, 2, 2)
    a5plt.hist1d(
        [xdata1, xdata2], xbins=10, legend=["1", "2"], weights=[wdata1, wdata2],
        xscale=scale, yscale=scale, xlabel="x-label", ylabel="y-label", axes=ax,
        )
    return fig

@pytest.mark.mpl_image_compare
@pytest.mark.parametrize(
    "scale",
    ["linear_pos", "linear_signed", "log_pos", "log_signed",],
    ids=["linear-pos", "linear-signed", "log-pos", "log-signed"],
    )
def test_hist2d(scale):
    xbins = np.sort(generate_data(10, mode=scale, seed=0))
    ybins = np.sort(generate_data(20, mode=scale, seed=1))
    xdata = generate_data(1000, mode=scale, seed=2)
    ydata = generate_data(1000, mode=scale, seed=3)
    wdata = generate_data(1000, mode=scale, seed=4)
    scale = "linear" if "linear" in scale else "log"

    fig, ax = a5plt.subplots()
    a5plt.hist2d(
        xdata, ydata, xbins=xbins, ybins=ybins, weights=wdata, xscale=scale,
        yscale=scale, cscale=scale, xlabel="x-label", ylabel="y-label",
        clabel="c-label", axes=ax,
        )
    return fig

@pytest.mark.mpl_image_compare
@pytest.mark.parametrize(
    "scale",
    ["linear_pos", "linear_signed", "log_pos", "log_signed",],
    ids=["linear-pos", "linear-signed", "log-pos", "log-signed"],
    )
def test_contour2d(scale):
    nx, ny = 200, 100
    x = np.sort(generate_data(nx, mode=scale, seed=0))
    y = np.sort(generate_data(ny, mode=scale, seed=1))
    x, y = np.meshgrid(x, y, indexing="ij")
    c = np.reshape(
        generate_data(nx*ny, mode=scale, seed=2), (nx, ny),
        )
    if "linear" in scale:
        c = np.cos(3*x) * np.sin(6*y)
        contours = [-0.8, -0.4, -0.2, 0.0, 0.2, 0.4, 0.8]
        if "pos" in scale:
            c += 1.0
            contours = [0.2, 0.4, 0.8, 1.0, 1.2, 1.4, 1.8]
    else:
        c = 1e7*np.exp(-5*np.log10((x+10)**2 + y**2))
        contours = np.array([1e1, 1e2, 1e3])
        if "signed" in scale:
            c -= 1e7*np.exp(-5*np.log10((x-10)**2 + y**2))
            contours = np.array([-1e3, -1e2, -1e1, 0, 1e1, 1e2, 1e3])
    scale = "linear" if "linear" in scale else "log"

    fig = a5plt.figure()
    ax = fig.add_subplot(1, 3, 1)
    a5plt.contour2d(
        x, y, c, contours, linestyles="solid", xscale=scale, yscale=scale,
        cscale=scale, xlabel="x-label", ylabel="y-label", clabel="c-label",
        labels=False, axes=ax,
        )
    ax = fig.add_subplot(1, 3, 2)
    a5plt.contour2d(
        x, y, c, contours, fill=True, extend_max="red", extend_min="blue",
        linestyles="solid", xscale=scale, yscale=scale, cscale=scale,
        xlabel="x-label", ylabel="y-label", clabel="c-label", axes=ax,
        )
    ax = fig.add_subplot(1, 3, 3)
    contours = 5
    a5plt.contour2d(
        x, y, c, contours, fill=True, linestyles="solid", xscale=scale,
        yscale=scale, cscale=scale, xlabel="x-label", ylabel="y-label",
        clabel="c-label", labels=["a", "b", "c", "d", "e"], axes=ax,
        )
    return fig

@pytest.mark.mpl_image_compare
def test_arrow2d():
    n = 10
    x1, y1, dx, dy, c = (
        generate_data(n, mode="linear_signed", seed=0),
        generate_data(n, mode="linear_signed", seed=1),
        generate_data(n, mode="linear_signed", seed=2),
        generate_data(n, mode="linear_signed", seed=3),
        generate_data(n, mode="linear_signed", seed=4),
        )

    fig, ax = a5plt.subplots()
    a5plt.arrows2d(
        x1[0]-0.5, y1[0], dx[0], dy[0], xlabel="x-label", ylabel="y-label",
        axes=ax,
        )
    a5plt.arrows2d(
        x1[0]+0.5, y1[0], dx[0], dy[0], "red", xlabel="x-label",
        ylabel="y-label", axes=ax,
        )
    a5plt.arrows2d(
        x1, y1, dx, dy, c, xlabel="x-label", ylabel="y-label", clabel="c-label",
        equal=True, axes=ax,
        )
    return fig
