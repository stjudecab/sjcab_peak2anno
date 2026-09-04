

no tests ran in 0.00s
ERROR: file or directory not found: tests



==================================== ERRORS ====================================
__________________ ERROR collecting test_peak2anno_package.py __________________
/home/bxu2/.local/lib/python3.6/site-packages/_pytest/python.py:599: in _importtestmodule
    mod = import_path(self.path, mode=importmode, root=self.config.rootpath)
/home/bxu2/.local/lib/python3.6/site-packages/_pytest/pathlib.py:533: in import_path
    importlib.import_module(module_name)
/usr/lib64/python3.6/importlib/__init__.py:126: in import_module
    return _bootstrap._gcd_import(name[level:], package, level)
<frozen importlib._bootstrap>:994: in _gcd_import
    ???
<frozen importlib._bootstrap>:971: in _find_and_load
    ???
<frozen importlib._bootstrap>:955: in _find_and_load_unlocked
    ???
<frozen importlib._bootstrap>:665: in _load_unlocked
    ???
/home/bxu2/.local/lib/python3.6/site-packages/_pytest/assertion/rewrite.py:162: in exec_module
    source_stat, co = _rewrite_test(fn, self.config)
/home/bxu2/.local/lib/python3.6/site-packages/_pytest/assertion/rewrite.py:366: in _rewrite_test
    co = compile(tree, strfn, "exec", dont_inherit=True)
E     File "/research/rgs01/home/clusterHome/bxu2/tgz/pypi/sjcab_peak2anno/test_peak2anno_package.py", line 3
E       from __future__ import annotations
E                                        ^
E   SyntaxError: future feature annotations is not defined
=========================== short test summary info ============================
ERROR test_peak2anno_package.py
!!!!!!!!!!!!!!!!!!!! Interrupted: 1 error during collection !!!!!!!!!!!!!!!!!!!!
1 error in 0.23s
Using CPython 3.13.3
Creating virtual environment at: .venv
error: Failed to fetch: `https://pypi.org/simple/matplotlib/`
  Caused by: Could not connect, are you offline?
  Caused by: Request failed after 3 retries
  Caused by: error sending request for url (https://pypi.org/simple/matplotlib/)
  Caused by: client error (Connect)
  Caused by: dns error: failed to lookup address information: Name or service not known
  Caused by: failed to lookup address information: Name or service not known

....                                                                     [100%]
4 passed in 0.08s

....                                                                     [100%]
4 passed in 0.02s

....                                                                     [100%]
4 passed in 0.02s

....                                                                     [100%]
4 passed in 0.03s
