from setuptools import setup

setup(
    name="WISP",
    version="3.0",
    description="Wenger Interferometry Software Package",
    author="Trey V. Wenger",
    author_email="tvwenger@gmail.com",
    packages=["wisp"],
    package_data={"wisp": ["logging.conf"]},
    install_requires=["astropy"],
)
