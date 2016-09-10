Cap
===

Cap is an open-source software library for modeling energy storage devices.
To get started, [checkout the documentation](https://cap.readthedocs.org).


Continuous Integration Status
-----------------------------

Provider      | Service            | Status
------------- | ------------------ | ------
Read the Docs | documentation      | [![Documentation Status](https://readthedocs.org/projects/cap/badge/?version=latest)](https://readthedocs.org/projects/cap/?badge=latest)
Travis CI     | unit tests         | [![Build Status](https://travis-ci.org/ORNL-CEES/Cap.svg?branch=master)](https://travis-ci.org/ORNL-CEES/Cap)
Codecov       | coverage           | [![codecov](https://codecov.io/gh/ORNL-CEES/Cap/branch/master/graph/badge.svg)](https://codecov.io/gh/ORNL-CEES/Cap)
Slack         | messaging          | [![Slack](https://img.shields.io/badge/Slack-%23cap-ff69b4.svg)](https://ornl-cees.slack.com/archives/cap)
Waffle        | project management | [![Stories in Ready](https://badge.waffle.io/ORNL-CEES/Cap.png?label=ready&title=Ready)](https://waffle.io/ORNL-CEES/Cap)


NEW - Escape dependency hell
----------------------------

No need to build Cap and its third-party libraries from source. Just pull the
latest image with Docker, run it, and launch numerical simulations via your
web browser. PyCap makes it possible to create and manipulate Cap
``EnergyStorageDevice`` objects and execute algorithms (``Charge``,
``Discharge``, etc.) from an interactive Python command-line interpreter,
without the need to constantly recompile.

[![Docker Pulls](https://img.shields.io/docker/pulls/dalg24/cap.svg)](https://hub.docker.com/r/dalg24/cap)
[![Image Layers](https://images.microbadger.com/badges/image/dalg24/cap.svg)](http://microbadger.com/images/dalg24/cap)


Authors
-------
* [Damien Lebrun-Grandie](https://github.com/dalg24)
* [Bruno Turcksin](https://github.com/rombur)
