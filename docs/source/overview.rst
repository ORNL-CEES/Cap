Getting started
===============

Overview
--------

`Cap <https://github.com/dalg24/cap>`_ is a library for modeling energy
storage devices.
Its core is implemented in C++ but Python wrappers are also available.

Cap provides:

1. energy storage device models
2. electrochemical measurement techniques

Guidelines for installation are :ref:`available <installation>`.
Note that it is **not** necessary to build Cap from source to use it.
For instructions on how to use Cap without installing it, refer to the
following :ref:`section <docker>`.


.. _docker:

Alternative to the full install procedure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
All you need is a working installation of Docker.
Follow the `Docker Engine installation guide
<https://docs.docker.com/engine/installation/>`_ for details on how to
install it on your machine.
It is supported on Linux, Cloud, Windows, and OS X.

The following command starts a Docker container with a Jupyter Notebook server
listening for HTTP connections on port 8888.
It mounts the present working directory, ``$PWD``, into the container at
``/notebooks``, which is set as the Jupyter Notebook startup folder.
It has pycap already installed on it and comes with a few notebooks as example.

::

    $ docker run --rm -it \
          -p 8888:8888 \
          -v $PWD:/notebooks \
          dalg24/cap

Open your web browser and follow ``http://localhost:8888``.
