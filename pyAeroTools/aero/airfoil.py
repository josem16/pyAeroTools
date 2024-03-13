import numpy as np
import matplotlib.pyplot as plt

def naca_four_digit_series(naca_number, chord=1.0, num_points=50):
    """
    Compute coordinates for airfoil upper and lower surfaces.

    Parameters
    ----------
    naca_number : int or str
        NACA four-digit series number.
    chord : float, optional, default=1.0
        Airfoil chord

    Returns
    -------
    airfoil : numpy array
        Airfoil upper and lower surface coordinates
        shape = (num_points, 4)
        columns = [x_upper, y_upper, x_lower, y_lower]

    """
    # Check input
    if isinstance(naca_number,int):
        naca_number = str(naca_number)
    msg = f'ERROR: Input naca number must be 4 digits. {len(naca_number)} digits given.'
    assert len(naca_number) == 4, msg

    # The first digit specifies the maximum camber (m) in percentage of the chord (airfoil length)
    # The second indicates the position of the maximum camber (p) in tenths of chord
    # The last two numbers provide the maximum thickness (t) of the airfoil in percentage of chord
    m = int(naca_number[0])/100.*chord
    p = int(naca_number[1])/10.*chord
    t = int(naca_number[2:])/100.*chord

    # Initialize x coordinates
    x_stations = np.linspace(0., chord, num=num_points)

    # Compute the mean camber line coordinates
    yc = []
    for x_station in x_stations:
        if x_station < p:
            yc_temp = m/p**2 * (2*p*x_station - x_station**2)
        else:
            yc_temp = m/(1-p)**2 * ((1-2*p) + 2*p*x_station - x_station**2)
        yc.append(yc_temp)

    # Compute the thickness distribution above and below the mean line
    yt = t/0.2 * (0.2969*np.sqrt(x_stations) - 0.1260*x_stations - 0.3516*x_stations**2 + 0.2843*x_stations**3 - 0.1015*x_stations**4)

    # Compute derivative of mean camber line w.r.t x_station
    dyc_dx = []
    for x_station in x_stations:
        if x_station < p:
            dyc_dx_temp = m/p**2 * (2*p - 2*x_station)
        else:
            dyc_dx_temp = m/(1-p)**2 * (2*p - 2*x_station)
        dyc_dx.append(dyc_dx_temp)

    # Compute airfoil upper and lower surface coordinates
    theta = np.arctan2(dyc_dx,1)
    x_upper = x_stations - yt*np.sin(theta)
    y_upper = yc + yt*np.cos(theta)
    x_lower = x_stations + yt*np.sin(theta)
    y_lower = yc - yt*np.cos(theta)
    airfoil = np.column_stack((x_upper, y_upper, x_lower, y_lower))
    return airfoil

def main():
    naca0012 = naca_four_digit_series('2317')
    fig=plt.figure()
    plt.plot(naca0012[:,0],naca0012[:,1])
    plt.plot(naca0012[:,2],naca0012[:,3])
    plt.axis('equal')
    plt.show()

if __name__ == '__main__':
    main()