# Unit tests

Standard-library unittest cases validate scientific edge cases and network-transfer logic through
mocks. They do not require external downloads. Source tests are tracked; temporary files are removed.
Run `python3 -m unittest discover -s tests/unit -v` from the root. Real tool accuracy is checked
separately by integration tests where available.
