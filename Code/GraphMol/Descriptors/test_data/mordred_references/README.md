# Mordred reference values

The `.json` files are what the osmordred tests read
(`Wrap/test_osmordred_ported_mordred.py`, `Wrap/test_osmordred_ic_mordred_flavor.py`),
so the tests need only the Python standard library.

The `.yaml` files are the original mordred reference data, kept for reference.
The tests do not read them. The JSON was generated from them and holds the same
values; comments above a YAML block are kept as a `"comment"` field in the JSON.
If you edit one format, update the other to match.
