from setuptools import setup
from setuptools.command.build import build
import subprocess
import shutil
import os
import re

class CustomBuild(build):
    # custom build commands
    def run(self):
        # run cmake
        try:
            subprocess.run(
                ["cmake","-B", "build_cmake"],
                check=True,
                capture_output=True)
        except subprocess.CalledProcessError as e:
            print( e.stdout)
            print(e.stderr)
            raise SystemExit(1)
        subprocess.run(["cmake", "--build", "build_cmake"])

        ## copy libira.so into site-packages/ira_mod/lib ::
        # libira.so path
        so_src = os.path.abspath(os.path.join('lib', 'libira.so'))
        # self.build_lib is location for site-packages
        build_dir_lib = os.path.join(self.build_lib, "ira_mod/lib")
        if not os.path.exists(build_dir_lib):
            os.makedirs(build_dir_lib)
        # destination path
        so_dst = os.path.join(build_dir_lib, 'libira.so')
        # copy
        shutil.copy2(so_src, so_dst)

        super().run()

# find the ira version string
with open(os.path.abspath("src/version.f90"), "r", encoding="utf-8") as f:
    for line in f:
        if( "string =" in line):
            vstr=line.strip()
# extract version format d.d.d
match = re.search( r"\d+\.\d+\.\d+", vstr)
ira_version=match.group()

setup(
    name="ira_mod",
    version=ira_version,
    description="",
    author="MAMMASMIAS Consortium",
    url="https://github.com/mammasmias/IterativeRotationsAssignments",
    license_expression="GPL-3.0-or-later OR Apache-2.0",
    packages=['ira_mod'],
    package_dir={'ira_mod': 'interface'},
    include_package_data=True,
    cmdclass={'build': CustomBuild},
    ext_modules=[],
)
