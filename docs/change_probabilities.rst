Change MC Probabilities
=======================

The probabilities of a batch of input files can be changed via

.. code-block:: python

    python3 /path/to/changeProbs.py -f /paths/to/file/dirs -i names

where :code:`/paths/to/file/dirs` are the paths to the parent directories
of the independent simulations, and :code:`names` are the names of the independent simulations.

.. note::
   it is recommended to add the :code:`-vv` argument to get extra verbosity to see what happens


For example, the invocation might be

.. code-block:: python

    python3 ~/MCFlow/mcflow/changeProbs.py -f *-0/*K/0.001MPa -i 1 2 3 4


To find out more options, ask for help


.. code-block:: python

    python3 /path/to/changeProbs.py --help
