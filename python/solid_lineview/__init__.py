# SPDX-License-Identifier: LGPL-3.0-or-later
"""solid_lineview: vector line drawings (SVG, PDF, EPS, PNG) of Geant4 geometry scenes.

The package is the renderer behind the ``solid_lineview`` command.  Typical use from Python::

    import solid_lineview as slv
    scene = slv.parse_prim('g4_0000.prim')            # a PrimScene
    slv.main(['g4_0000.prim', '--mode', 'section', '--section', 'x=0', '-o', 'side.svg'])

``python3 -m solid_lineview ...`` is the same as the ``solid_lineview`` command.  Everything lives
in :mod:`solid_lineview.core`; the names below are the supported interface, helpers with a leading
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
    main,
    parse_prim,
    render_section,
    write_outputs,
)
