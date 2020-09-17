from yaml import load, Loader


def get_velocity(cmission, mission_management):
    """
    Get velocity of a mission from MissionManagement.yaml
    """
    yaml = load(open(str(mission_management)), Loader=Loader)
    velocity = yaml[cmission]['VELOCITY']
    if velocity is not None:
        return velocity
    else:
        raise ValueError("velocity not found for mission %s in %s" % (cmission, mission_management))


def get_deltat(cmission, mission_management):
    """
    Get deltaT of a mission from file MissionManagement.yaml
    """
    yaml = load(open(str(mission_management)), Loader=Loader)
    deltat = yaml[cmission]['DELTA_T']
    if deltat is not None:
        return deltat
    else:
        raise ValueError("deltat not found for mission %s in %s" % (cmission, mission_management))
        