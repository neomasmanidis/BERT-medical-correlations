import os


if not os.path.exists('data dump'):
    print("Creating folder 'data dump'..")
    os.mkdir('data dump')

import Prep
import Curate
import GrabNames
