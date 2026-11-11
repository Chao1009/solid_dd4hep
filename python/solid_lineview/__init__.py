# SPDX-License-Identifier: LGPL-3.0-or-later
"""solid_lineview: vector line drawings (SVG, PDF, EPS, PNG) of Geant4 geometry scenes.

The package is the renderer behind the ``solid_lineview`` command.  Typical use from Python::

    import solid_lineview as slv
    scene = slv.load_scene('scene.npz')               # a PrimScene, from a scene.npz or a .prim file
    slv.main(['scene.npz', '--mode', 'section', '--section', 'x=0', '-o', 'side.svg'])

``python3 -m solid_lineview ...`` is the same as the ``solid_lineview`` command.  Everything lives
in :mod:`solid_lineview.core` (renderer) and :mod:`solid_lineview.scene` (scene.npz reader and
.prim writer); the names below are the supported interface, helpers with a leading
underscore are reachable as ``solid_lineview.core._name``.
"""
from .core import (  # noqa: F401
    __version__,
    Camera,
    DrawList,
    FieldMap,
    Mesh,
    PrimScene,
    Renderer,
    Solid,
    build_parser,
    cairo_draw_items,
    load_scene,
    main,
    parse_prim,
    render_section,
    svg_items,
    write_outputs,
)
