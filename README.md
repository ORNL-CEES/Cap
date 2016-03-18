Cap
===
[![GitHub license](https://img.shields.io/github/license/ORNL-CEES/Cap.svg)]()

Cap is an open-source software library for modeling energy storage devices.
To get started, [checkout the documentation](https://cap.readthedocs.org).

Continuous Integration Status
-----------------------------

Provider      | Service       | Status
------------- | ------------- | ------
Read the Docs | documentation | [![Documentation Status](https://readthedocs.org/projects/cap/badge/?version=latest)](https://readthedocs.org/projects/cap/?badge=latest)
Travis CI     | unit tests    | [![Build Status](https://travis-ci.org/ORNL-CEES/Cap.svg?branch=master)](https://travis-ci.org/ORNL-CEES/Cap)
Codecov       | coverage      | [![codecov.io](https://codecov.io/github/ORNL-CEES/Cap/coverage.svg?branch=master)](https://codecov.io/github/ORNL-CEES/Cap?branch=master)
Slack         | messaging     | [![Slack](https://img.shields.io/badge/Slack-%23cap-ff69b4.svg)](https://ornl-cees.slack.com/archives/cap)


NEW - Escape dependency hell
----------------------------

No need to build Cap and its third-party libraries from source. Just pull the
latest image with Docker, run it, and launch numerical simualtions via your
web browser. PyCap makes it possible to create and manipulate Cap
``EnergyStorageDevice`` objects and execute algorithms (``Charge``,
``Discharge``, etc.) from an interactive Python command-line interpreter,
without the need to constantly recompile.

[![Docker Pulls](https://img.shields.io/docker/pulls/dalg24/cap.svg)](https://hub.docker.com/r/dalg24/cap)
[![ImageLayers Size](https://img.shields.io/imagelayers/image-size/dalg24/cap-stack/latest.svg)]()

Authors
-------
* [Damien Lebrun-Grandie](https://github.com/dalg24)
* [Bruno Turcksin](https://github.com/rombur)
