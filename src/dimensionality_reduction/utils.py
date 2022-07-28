import pkg_resources

def check_version(pkg):
    """Get the version of a Python package that may or may not be installed."""
    try:
        return pkg_resources.get_distribution(pkg).version
    except pkg_resources.DistributionNotFound:
        return "ModuleNotFound"
