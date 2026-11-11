"""Shared fixtures: one built toy project per session, a fresh Frame per test."""
from __future__ import annotations

import pytest

from coords import Frame
from coords.build import build_all
from coords.config import load_config
from coords.testing import Toy, make_toy_project


@pytest.fixture(scope="session")
def toy(tmp_path_factory) -> Toy:
    t = make_toy_project(tmp_path_factory.mktemp("toy"))
    build_all(load_config(t.config), verbose=False)
    return t


@pytest.fixture()
def fr(toy) -> Frame:
    """A fresh Frame (empty hole log) over the built toy maps."""
    return Frame.load(toy.maps)
