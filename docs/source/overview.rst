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

    $ docker run -d \
          -p 8888:8888 \
          -v $PWD:/notebooks \
          dalg24/cap

Open your web browser and follow ``http://<ip_address>:8888``
where ``<ip_address>`` is the IP address of the machine with the Docker daemon
running.

- On OS X, the IP address can be obtained by:

::

    # assuming that you are running docker on the "default" VM
    $ docker-machine ip default
    192.168.99.100 # <- this is the machine ip address

In that case you would copy/paste the url http://192.168.99.100:8888 into
your browser address bar.

- On Linux, you may access the notebook server from the browser using
  http://localhost:8888. You may have use ``sudo`` to run the ``docker``
  command. Please refer to the Docker documentation on how to create a
  ``docker`` group and add your user to it if you want to avoid using
  ``sudo`` each time.
