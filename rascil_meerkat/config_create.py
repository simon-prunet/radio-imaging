# pylint: disable=invalid-name,too-many-arguments,too-many-locals

"""
Create configuration.
"""

import logging
from typing import Union

import numpy
import pkg_resources
from astropy import units
from astropy.coordinates import EarthLocation

from ska_sdp_datamodels.configuration.config_coordinate_support import (
    ecef_to_enu,
    lla_to_ecef,
)
from ska_sdp_datamodels.configuration.config_model import Configuration

log = logging.getLogger("data-models-logger")


ANTENNA_FILES = {
    "LOWBD2": pkg_resources.resource_stream(
        __name__, "example_antenna_files/ska/LOWBD2.csv"
    ).name,
    "LOWBD2-CORE": pkg_resources.resource_stream(
        __name__, "example_antenna_files/ska/LOWBD2-CORE.csv"
    ).name,
    "LOW": pkg_resources.resource_stream(
        __name__, "example_antenna_files/ska/LOW_SKA-TEL-SKO-0000422_Rev3.txt"
    ).name,
    "LOWR3": pkg_resources.resource_stream(
        __name__, "example_antenna_files/ska/LOW_SKA-TEL-SKO-0000422_Rev3.txt"
    ).name,
    "LOW-AA0.5": pkg_resources.resource_stream(
        __name__, "example_antenna_files/ska/SKA1-LOW-AA0.5-v1.1_AA0p5.txt"
    ).name,
    "MID": pkg_resources.resource_stream(
        __name__, "example_antenna_files/ska/ska1mid.cfg"
    ).name,
    "MIDR5": pkg_resources.resource_stream(
        __name__, "example_antenna_files/ska/ska1mid.cfg"
    ).name,
    "MID-AA0.5": pkg_resources.resource_stream(
        __name__, "example_antenna_files/ska/ska1mid.cfg"
    ).name,
    "MEERKAT": pkg_resources.resource_stream(
        __name__, "example_antenna_files/other/mkat.cfg"
    ).name,
    "MEERKAT+": pkg_resources.resource_stream(
        __name__, "example_antenna_files/other/mkatplus.cfg"
    ).name,
    "ASKAP": pkg_resources.resource_stream(
        __name__, "example_antenna_files/other/askap.cfg"
    ).name,
    "LOFAR": pkg_resources.resource_stream(
        __name__, "example_antenna_files/other/lofar.cfg"
    ).name,
    "VLAA": pkg_resources.resource_stream(
        __name__, "example_antenna_files/other/vlaa_local.csv"
    ).name,
    "VLAA_north": pkg_resources.resource_stream(
        __name__, "example_antenna_files/other/vlaa_local.csv"
    ).name,
}


def _limit_rmax(antxyz, diameters, names, mounts, rmax):
    """
    Select antennas with radius from centre < rmax

    :param antxyz: Geocentric coordinates
    :param diameters: diameters in metres
    :param names: Names
    :param mounts: Mount types
    :param rmax: Maximum radius (m)
    :return: sub-selected antxyz, diameters, names, mounts
    """
    if rmax is not None:
        lantxyz = antxyz - numpy.average(antxyz, axis=0)
        r = numpy.sqrt(
            lantxyz[:, 0] ** 2 + lantxyz[:, 1] ** 2 + lantxyz[:, 2] ** 2
        )

        antxyz = antxyz[r < rmax]
        log.debug(
            "create_configuration_from_file: Maximum radius "
            "%.1f m includes %d antennas/stations",
            rmax,
            antxyz.shape[0],
        )
        diameters = diameters[r < rmax]
        names = numpy.array(names)[r < rmax]
        mounts = numpy.array(mounts)[r < rmax]
    else:
        log.debug(
            "create_configuration_from_file: %d antennas/stations",
            antxyz.shape[0],
        )
    return antxyz, diameters, names, mounts


def _find_vptype_from_name(names, match: Union[str, dict] = "unknown"):
    """
    Determine voltage pattern type from name using a dictionary

    There are two modes:

    If match is a dict, then the antenna/station names are matched.
    An example of match would be: d={"M0":"MEERKAT", "SKA":"MID"}
    The test if whether the key e.g. M0 is in the
    antenna/station name e.g. M053

    If match is a str then the returned array is filled with that value.

    :param names: list/array of antenna names
    :param match: dict or str of voltage patter to be matched with names
    :return: voltage pattern type
    """
    if isinstance(match, dict):
        vp_types = numpy.repeat("unknown", len(names))
        for item in match:
            for i, name in enumerate(names):
                if item in name:
                    vp_types[i] = match.get(item)
    elif isinstance(match, str):
        vp_types = numpy.repeat(match, len(names))
    else:
        raise ValueError("match must be str or dict")

    return vp_types


def create_configuration_from_file(
    antfile: str,
    location: EarthLocation = None,
    mount: str = "azel",
    names: str = "%d",
    vp_type: Union[str, dict] = "Unknown",
    diameter=35.0,
    rmax=None,
    name="",
    skip=1,
    ecef=True,
) -> Configuration:
    """
    Define configuration from a text file

    :param antfile: Antenna file name
    :param location: Earthlocation of array
    :param mount: mount type: 'azel', 'xy', 'equatorial'
    :param names: Antenna names e.g. "VLA%d"
    :param vp_type: String or rule to map name to voltage pattern type
    :param diameter: Effective diameter of station or antenna
    :param rmax: Maximum distance from array centre (m)
    :param name: Name of array
    :param skip: Antennas/stations to skip (int)
    :param ecef: Configuration file format: ECEF or local-xyz;
                 if ECEF, it will be converted to ENU coordinates
    :return: Configuration
    """
    antxyz = numpy.genfromtxt(antfile, delimiter=",")
    if not antxyz.shape[1] == 3:
        raise ValueError(f"Antenna array has wrong shape {antxyz.shape}")

    nants = antxyz.shape[0]
    if ecef:
        antxyz = ecef_to_enu(location, antxyz)

    diameters = diameter * numpy.ones(nants)
    anames = [names % ant for ant in range(nants)]
    mounts = numpy.repeat(mount, nants)
    antxyz, diameters, anames, mounts = _limit_rmax(
        antxyz, diameters, anames, mounts, rmax
    )

    antxyz = antxyz[::skip]
    diameters = diameters[::skip]
    anames = anames[::skip]
    mounts = mounts[::skip]

    fc = Configuration.constructor(
        location=location,
        names=anames,
        mount=mounts,
        xyz=antxyz,
        vp_type=_find_vptype_from_name(anames, vp_type),
        diameter=diameters,
        name=name,
    )
    return fc


def create_configuration_from_MIDfile(
    antfile: str,
    location=None,
    mount: str = "azel",
    vp_type: Union[str, dict] = "Unknown",
    rmax=None,
    name="",
    skip=1,
    ecef=True,
) -> Configuration:
    """
    Define configuration from a SKA MID format file

    :param antfile: Antenna file name
    :param location: EarthLocation of array
    :param mount: mount type: 'azel', 'xy'
    :param vp_type: String or rule to map name to voltage pattern type
    :param rmax: Maximum distance from array centre (m)
    :param name: Name of array
    :param skip: Antennas/stations to skip
    :param ecef: Configuration file format: ECEF or local-xyz;
                 if ECEF, it will be converted to ENU coordinates
    :return: Configuration
    """

    antxyz = numpy.genfromtxt(
        antfile, skip_header=5, usecols=[0, 1, 2], delimiter=" "
    )

    nants = antxyz.shape[0]
    if not antxyz.shape[1] == 3:
        raise ValueError(f"Antenna array has wrong shape {antxyz.shape}")
    if ecef:
        antxyz = ecef_to_enu(location, antxyz)

    anames = numpy.genfromtxt(
        antfile, dtype="str", skip_header=5, usecols=[4], delimiter=" "
    )
    mounts = numpy.repeat(mount, nants)
    diameters = numpy.genfromtxt(
        antfile, dtype="str", skip_header=5, usecols=[3], delimiter=" "
    ).astype("float")

    antxyz, diameters, anames, mounts = _limit_rmax(
        antxyz, diameters, anames, mounts, rmax
    )

    antxyz = antxyz[::skip]
    diameters = diameters[::skip]
    anames = anames[::skip]
    mounts = mounts[::skip]

    fc = Configuration.constructor(
        location=location,
        names=anames,
        mount=mounts,
        xyz=antxyz,
        vp_type=_find_vptype_from_name(anames, vp_type),
        diameter=diameters,
        name=name,
    )

    return fc


def create_configuration_from_LLAfile(
    antfile: str,
    location: EarthLocation = None,
    mount: str = "azel",
    names: str = "%d",
    vp_type: Union[str, dict] = "Unknown",
    diameter=35.0,
    alt=0.0,
    rmax=None,
    name="",
    skip=1,
    ecef=False,
) -> Configuration:
    """
    Define configuration from a longitude-latitude file

    :param antfile: Antenna file name
    :param location: Earthlocation of array
    :param mount: mount type: 'azel', 'xy', 'equatorial'
    :param names: Antenna names e.g. "VLA%d"
    :param vp_type: string or rule to map name to voltage pattern type
    :param diameter: Effective diameter of station or antenna
    :param alt: The altitude assumed
    :param rmax: Maximum distance from array centre (m)
    :param name: Name of array
    :param skip: Antennas/stations to skip
    :param ecef: Configuration file format: ECEF or local-xyz
                 if ECEF, it will be converted to ENU coordinates
    :return: Configuration
    """

    antxyz = numpy.genfromtxt(antfile, delimiter=",")

    nants = antxyz.shape[0]

    lon, lat = antxyz[:, 0], antxyz[:, 1]
    x, y, z = lla_to_ecef(lat * units.deg, lon * units.deg, alt)
    antxyz = numpy.stack((x, y, z), axis=1)

    if ecef:
        antxyz = ecef_to_enu(location, antxyz)

    diameters = diameter * numpy.ones(nants)
    anames = [names % ant for ant in range(nants)]
    mounts = numpy.repeat(mount, nants)
    antxyz, diameters, anames, mounts = _limit_rmax(
        antxyz, diameters, anames, mounts, rmax
    )

    antxyz = antxyz[::skip]
    diameters = diameters[::skip]
    anames = anames[::skip]
    mounts = mounts[::skip]

    fc = Configuration.constructor(
        location=location,
        names=anames,
        mount=mounts,
        xyz=antxyz,
        vp_type=_find_vptype_from_name(anames, vp_type),
        diameter=diameters,
        name=name,
    )
    return fc


# pylint: disable=too-many-branches,too-many-statements
def create_named_configuration(
    name: str = "LOWBD2", **kwargs
) -> Configuration:
    """
    Create standard configurations e.g. LOWBD2, MIDBD2

    Possible configurations are::
        LOWBD2
        LOWBD2-CORE
        LOW == LOWR3
        LOW-AA0.5
        MID == MIDR5
        MID-AA0.5
        MEERKAT
        MEERKAT+
        ASKAP
        LOFAR
        VLAA
        VLAA_north

    :param name: name of Configuration e.g. MID, LOW, LOFAR, VLAA, ASKAP
    :param rmax: Maximum distance of station from the average (m)
    :return: Configuration object

    For LOWBD2, setting rmax gives the following number of stations
    100.0       13
    300.0       94
    1000.0      251
    3000.0      314
    10000.0     398
    30000.0     476
    100000.0    512
    """
    low_location = EarthLocation(
        lon=116.76444824 * units.deg,
        lat=-26.824722084 * units.deg,
        height=300.0,
    )
    mid_location = EarthLocation(
        lon=21.443803 * units.deg,
        lat=-30.712925 * units.deg,
        height=1053.000000,
    )
    meerkat_location = EarthLocation(
        lon=21.44388889 * units.deg, lat=-30.7110565 * units.deg, height=1086.6
    )
    if name == "LOWBD2":
        location = low_location
        log.debug(
            "create_named_configuration: %s\n\t%s\n\t%s",
            name,
            location.geocentric,
            location.geodetic,
        )
        fc = create_configuration_from_file(
            antfile=ANTENNA_FILES[name],
            location=location,
            mount="XY",
            names="LOWBD2_%d",
            vp_type="LOW",
            diameter=35.0,
            name=name,
            ecef=False,
            **kwargs,
        )
    elif name == "LOWBD2-CORE":
        location = low_location
        log.debug(
            "create_named_configuration: %s\n\t%s\n\t%s",
            name,
            location.geocentric,
            location.geodetic,
        )
        fc = create_configuration_from_file(
            antfile=ANTENNA_FILES[name],
            vp_type="LOW",
            location=location,
            mount="XY",
            names="LOWBD2_%d",
            diameter=35.0,
            name=name,
            ecef=False,
            **kwargs,
        )
    elif name in ["LOW", "LOWR3"]:
        location = low_location
        log.debug(
            "create_named_configuration: %s\n\t%s\n\t%s",
            name,
            location.geocentric,
            location.geodetic,
        )
        fc = create_configuration_from_LLAfile(
            antfile=ANTENNA_FILES[name],
            location=location,
            mount="XY",
            names="SKA%03d",
            vp_type="LOW",
            diameter=38.0,
            alt=300.0,
            name=name,
            ecef=True,
            **kwargs,
        )
    elif name == "LOW-AA0.5":
        location = EarthLocation(
            lon=116.69345390 * units.deg,
            lat=-26.86371635 * units.deg,
            height=300.0,
        )

        log.debug(
            "create_named_configuration: %s\n\t%s\n\t%s",
            name,
            location.geocentric,
            location.geodetic,
        )
        fc = create_configuration_from_LLAfile(
            antfile=ANTENNA_FILES[name],
            location=location,
            mount="XY",
            names="LOW_AA0.5_%d",
            vp_type="LOW",
            diameter=38.0,
            alt=300.0,
            name=name,
            ecef=True,
            **kwargs,
        )
        # pylint: disable=duplicate-code
        names = [
            "S008‐1",
            "S008‐2",
            "S009‐1",
            "S009‐2",
            "S010‐1",
            "S010‐2",
        ]
        fc["names"].data = names

    elif name in ["MID", "MIDR5"]:
        location = mid_location
        log.debug(
            "create_named_configuration: %s\n\t%s\n\t%s",
            name,
            location.geocentric,
            location.geodetic,
        )
        fc = create_configuration_from_MIDfile(
            antfile=ANTENNA_FILES[name],
            vp_type={"M0": "MEERKAT", "SKA": "MID"},
            mount="azel",
            name=name,
            location=location,
            **kwargs,
        )
    elif name == "MID-AA0.5":
        location = mid_location
        log.debug(
            "create_named_configuration: %s\n\t%s\n\t%s",
            name,
            location.geocentric,
            location.geodetic,
        )
        fc = create_configuration_from_MIDfile(
            antfile=ANTENNA_FILES[name],
            vp_type={"M0": "MEERKAT", "SKA": "MID"},
            mount="azel",
            name=name,
            location=location,
            **kwargs,
        )
        # See Ben Mort's comment on
        # https://confluence.skatelescope.org/display/SE/Requirements+for+RCAL+pipelines+for+MID+and+LOW+AA0.5
        names = ["SKA001", "SKA100", "SKA036", "SKA063"]
        fc = select_configuration(fc, names)
    elif name == "MEERKAT":
        location = meerkat_location
        log.debug(
            "create_named_configuration: %s\n\t%s\n\t%s",
            name,
            location.geocentric,
            location.geodetic,
        )
        fc = create_configuration_from_MIDfile(
            antfile=ANTENNA_FILES[name],
            vp_type={"m0": "MEERKAT"},
            mount="ALT-AZ",
            name=name,
            location=location,
            **kwargs,
        )
    elif name == "MEERKAT+":
        location = meerkat_location
        log.debug(
            "create_named_configuration: %s\n\t%s\n\t%s",
            name,
            location.geocentric,
            location.geodetic,
        )
        fc = create_configuration_from_MIDfile(
            antfile=ANTENNA_FILES[name],
            vp_type={"m0": "MEERKAT", "s0": "MID"},
            mount="ALT-AZ",
            name=name,
            location=location,
            **kwargs,
        )
    elif name == "ASKAP":
        location = EarthLocation(
            lon=+116.6356824 * units.deg,
            lat=-26.7013006 * units.deg,
            height=377.0 * units.m,
        )
        log.debug(
            "create_named_configuration: %s\n\t%s\n\t%s",
            name,
            location.geocentric,
            location.geodetic,
        )
        fc = create_configuration_from_MIDfile(
            antfile=ANTENNA_FILES[name],
            vp_type="ASKAP",
            mount="equatorial",
            name=name,
            location=location,
            **kwargs,
        )
    elif name == "LOFAR":
        location = EarthLocation(
            x=3826923.9 * units.m, y=460915.1 * units.m, z=5064643.2 * units.m
        )
        log.debug(
            "create_named_configuration: %s\n\t%s\n\t%s",
            name,
            location.geocentric,
            location.geodetic,
        )

        if not kwargs.get("meta", False) is False:
            raise ValueError("For LOFAR, do not use the 'meta' input kwarg.")

        fc = create_configuration_from_MIDfile(
            antfile=ANTENNA_FILES[name],
            location=location,
            mount="XY",
            vp_type="LOFAR",
            name=name,
            **kwargs,
        )
    elif name == "VLAA":
        location = EarthLocation(
            lon=-107.6184 * units.deg, lat=34.0784 * units.deg, height=2124.0
        )
        log.debug(
            "create_named_configuration: %s\n\t%s\n\t%s",
            name,
            location.geocentric,
            location.geodetic,
        )
        fc = create_configuration_from_file(
            antfile=ANTENNA_FILES[name],
            location=location,
            mount="AZEL",
            names="VLA_%d",
            vp_type="VLA",
            diameter=25.0,
            name=name,
            ecef=False,
            **kwargs,
        )
    elif name == "VLAA_north":
        location = EarthLocation(
            lon=-107.6184 * units.deg, lat=90.000 * units.deg, height=0.0
        )
        log.debug(
            "create_named_configuration: %s\n\t%s\n\t%s",
            name,
            location.geocentric,
            location.geodetic,
        )
        fc = create_configuration_from_file(
            antfile=ANTENNA_FILES[name],
            location=location,
            mount="AZEL",
            names="VLA_%d",
            vp_type="VLA",
            diameter=25.0,
            name=name,
            ecef=False,
            **kwargs,
        )
    else:
        raise ValueError(f"No such Configuration {name}")
    return fc


def select_configuration(config, names=None):
    """
    Select a subset of antennas based on their "names"

    :param config: Configuration object
    :param names: names of antennas to sub-select
    :return: new Configuration with selected antennas
    """

    if names is None:
        return config

    names = numpy.array(names)
    ind = []
    for iname, name in enumerate(config.names.data):
        for aname in names:
            if aname.strip() == name.strip():
                ind.append(iname)

    if len(ind) == 0:
        raise ValueError(f"No antennas selected using names {names}")

    fc = Configuration.constructor(
        location=config.location,
        names=config.names[ind],
        mount=config.mount[ind],
        xyz=config.xyz[ind],
        vp_type=config.vp_type[ind],
        diameter=config.diameter[ind],
        name=config.name,
    )
    return fc


def decimate_configuration(config, start=0, stop=None, skip=1):
    """
    Decimate a configuration

    :param config: Configuration object
    :param start: Start of selection (int, index)
    :param stop: End of selection (int, index);
                 if None, include very last element
    :param skip: Antennas/stations to skip (int)
                 e.g if skip=2, it keeps every second antenna
                 counting from start+1
                 if skip=3, it keeps every third antenna
                 counting from start+1
    :return: decimated Configuration (new object)
    """
    fc = Configuration.constructor(
        location=config.location,
        names=config.names[start:stop:skip],
        mount=config.mount[start:stop:skip],
        xyz=config.xyz[start:stop:skip],
        vp_type=config.vp_type[start:stop:skip],
        diameter=config.diameter[start:stop:skip],
        name=config.name,
    )
    return fc
