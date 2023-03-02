from os.path import join
from library.TEST_UTILS import getfilelist

for file in sorted(getfilelist("./test")):
    file_path = join('./test', file)
    with open(file_path) as f:
        code_to_run = f.read()
    exec(code_to_run)
