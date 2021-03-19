import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize as opt
import math_eq

# class made to seperate rows of data in data file
class CepheidStar:
    def __init__(self, row):
        self.star_name = row[0]
        self.parallax = float(row[1])
        self.p_err = float(row[2])
        self.period = float(row[3])
        self.v_band = float(row[4])
        self.i_band = float(row[5])
        self.v_extinction = float(row[6])
        self.av_err = float(row[7])

        # equations to calculate various aspects to plot PL relationship graph
        # see math_eq for more detailed explanation
        # actual variables/data I want is listed in same order from class
        self.distance_in_parsecs = math_eq.get_distance_in_parsecs(self.parallax)

        self.M_v = math_eq.get_absolute_magnitude_visible(self.v_band,
                                                          self.distance_in_parsecs,
                                                          self.v_extinction, )

        self.M_i = math_eq.get_absolute_magnitude_infra_red(self.i_band,
                                                            self.distance_in_parsecs,
                                                            self.v_extinction)

        self.log_period = math_eq.get_log_period(self.period)

        self.distance_err = math_eq.get_distance_err(self.p_err, self.parallax)

        self.err = math_eq.get_absolute_magnitude_visible_error(self.distance_err,
                                                                self.distance_in_parsecs,
                                                                self.av_err)

        self.err_i = math_eq.get_absolute_magnitude_infra_red_error(self.distance_err,
                                                                    self.distance_in_parsecs,
                                                                    self.av_err)
# class used to seperate rows of data
class Galaxy:
    def __init__(self, galaxy_row, star_data):
        self.galaxy_name = galaxy_row[0]
        self.velocity = float(galaxy_row[1])
        self.A_VMW = float(galaxy_row[2])
        # making an array of stars from the data hstgalaxies in star_data
        self.stars = [HstGalaxyCepheids(star) for star in star_data]


# class used to seperate rows of data
class HstGalaxyCepheids:
    def __init__(self, row):
        self.cepheid_name = row[0]
        self.log_P = float(row[1])
        self.m_V = float(row[2])
        self.m_IR = float(row[3])

def get_Hubbles_constant():
    # load cephieds data as strings
    mw_data = np.loadtxt("data/MW_Cepheids.dat", dtype=str)

    # create star list
    cepheid_stars = []

    # create star instances
    for row in mw_data:
        current_star = CepheidStar(row)
        cepheid_stars.append(current_star)

    # Set graph labels
    plt.xlabel("Log Period")
    plt.ylabel("Absolute Magnitude")
    plt.title("PL relation")

    # Plot star points one at a time
    for star in cepheid_stars:
        # infra red plot
        plt.plot(star.log_period,
                 star.M_i,
                 color='red',
                 marker='x',
                 markersize=5,
                 linestyle="None")

        # infra red error
        plt.errorbar(star.log_period,
                     star.M_i,
                     yerr=star.err_i,
                     color='red',
                     linewidth=2,
                     linestyle="None")

        # visible
        plt.plot(star.log_period,
                 star.M_v,
                 color='green',
                 marker='x',
                 markersize=5,
                 linestyle="None")

        # visible error
        plt.errorbar(star.log_period,
                     star.M_v,
                     yerr=star.err,
                     color='green',
                     linewidth=2,
                     linestyle="None")

    # making data and errors into lists with multiple items
    xdata = [star.log_period for star in cepheid_stars]
    visible_y_data = [star.M_v for star in cepheid_stars]
    ir_y_data = [star.M_i for star in cepheid_stars]
    err_i = [star.err_i for star in cepheid_stars]
    err_v = [star.err for star in cepheid_stars]

    # using curve fit to find fit and intercept, and covariance to find error
    alpha_beta_i, covarianceI = opt.curve_fit(math_eq.get_period_luminosity_relation,
                                              xdata,
                                              ir_y_data,
                                              sigma=err_i)
    alpha_beta_v, covarianceV = opt.curve_fit(math_eq.get_period_luminosity_relation,
                                              xdata,
                                              visible_y_data,
                                              sigma=err_v)
    # setting alpha and beta into individual varibales from point in a list
    ir_alpha = alpha_beta_i[0]
    visible_alpha = alpha_beta_v[0]
    ir_beta = alpha_beta_i[1]
    visible_beta = alpha_beta_v[1]

    # calculating errors, square rooting and making into a diagonal matrix for covariance
    errorsI = np.sqrt(np.diag(covarianceI))
    errorsV = np.sqrt(np.diag(covarianceV))
    # making erros specific variables
    errorI_alpha = errorsI[0]
    errorV_alpha = errorsV[0]
    errorI_beta = errorsI[1]
    errorV_beta = errorsV[1]

    print("Alpha for IR is ", round(ir_alpha, 3), " ± ", round(errorI_alpha, 3),
          "and beta is ", round(ir_beta, 3), " ± ", round(errorI_beta, 3))
    print("Alpha for visible is ", round(visible_alpha, 3), " ± ",
          round(errorV_alpha, 3), "and beta is ", round(visible_beta, 3), " ± ",
          round(errorV_beta, 3))
    # return row of 100 evenly spaces data of the minimum and max od the xdata/log period
    x = np.linspace(min(xdata),
                    max(xdata),
                    100)
    # pltting line of best fit
    plt.plot(x,
             x * alpha_beta_i[0] + alpha_beta_i[1],
             color='red',
             label="Infrared")

    plt.plot(x,
             x * alpha_beta_v[0] + alpha_beta_v[1],
             color='green',
             label="Visible")

    chiIR = math_eq.chisq(ir_y_data, np.array(xdata) * alpha_beta_i[0]
                          + alpha_beta_i[1], np.array(err_i))
    # calculating chi squared using math_eq
    chiVR = math_eq.chisq(visible_y_data, np.array(xdata) * alpha_beta_v[0]
                          + alpha_beta_v[1], np.array(err_v))

    print("chisq for IR is ", round(chiIR, 3))
    print("chisq for visible is ", round(chiVR, 3))
    # calculating reduced chi sqaured
    print("Reduced infrared Chi Square is", round(chiIR / (len(xdata) - 2), 3))
    print("Reduced visible Chi Square is", round(chiVR / (len(xdata) - 2), 3))

    # show the generated graph and legend
    plt.legend()
    plt.show()

    # Loading the data into classes
    galaxy_data = np.loadtxt("data/galaxy_data.dat", dtype=str)

    galaxy_list = []
    length = []
    
    # for loop that counts through the number of galaxies in glalaxy data
    # loads star_data running through hst_gal of any number .
    # the current galaxy it looks through is a particualr galaxy with particualar stars.
    # both length and galaxy list are added to for the number of stars
    # to make an average and to calculate distance.
    for i, galaxy in enumerate(galaxy_data):
        file_name = "data/hst_gal{}_cepheids.dat".format(i + 1)
        star_data = np.loadtxt(file_name, dtype=str)
        current_galaxy = Galaxy(galaxy, star_data)
        length.append(len(star_data))
        galaxy_list.append(current_galaxy)

    # labelling graph

    plt.xlabel("Distance (Mpc)")
    plt.ylabel("Recessional Velocity (km/s)")
    plt.title("Recessional velocity against distance of all stars")

    y_data_list = []
    galaxy_velocity_list = []
    xV_data_list = []
    xIR_data_list = []
    xV_data_error_list = []
    xIR_data_error_list = []
    
    # going through galaxies in galaxy_list
    for galaxy in galaxy_list:
        # making a list of velocities
        velocity = [galaxy.velocity for galaxy in galaxy_list]
        # going through stars in each galaxy
        for star in galaxy.stars:
            # calculating distances
            distance_ir = math_eq.get_distance_ir(star.m_IR,
                                                  ir_alpha,
                                                  star.log_P,
                                                  ir_beta,
                                                  galaxy.A_VMW)
            distance_visible = math_eq.get_distance_visible(star.m_V,
                                                            visible_alpha,
                                                            star.log_P,
                                                            visible_beta,
                                                            galaxy.A_VMW)
            # calcualting errors
            error_absolute_mag_IR = math_eq.error_propogation(np.array(distance_ir),
                                                              errorI_alpha,
                                                              errorI_beta)
            error_absolute_mag_V = math_eq.error_propogation(np.array(distance_visible),
                                                             errorV_alpha,
                                                             errorV_beta)

            # making data lists by adding to them twice (theres the same
            # velocities for IR and V)
            galaxy_velocity_list.append(galaxy.velocity)
            y_data_list.append(galaxy.velocity)
            y_data_list.append(galaxy.velocity)
            # adding to the other lists to plot
            xIR_data_list.append(distance_ir)
            xV_data_list.append(distance_visible)
            xV_data_error_list.append(error_absolute_mag_V)
            xIR_data_error_list.append(error_absolute_mag_IR)

    # Plot points
    plt.plot(xIR_data_list,
             galaxy_velocity_list,
             color='red',
             marker='x',
             markersize=5,
             linestyle="None")

    # plot error
    plt.errorbar(xIR_data_list,
                 galaxy_velocity_list,
                 xerr=xIR_data_error_list,
                 color='red',
                 label="Infrared",
                 linewidth=2,
                 linestyle="None")

    plt.plot(xV_data_list,
             galaxy_velocity_list,
             color='green',
             marker='x',
             markersize=5,
             linestyle="None")

    plt.errorbar(xV_data_list,
                 galaxy_velocity_list,
                 xerr=xV_data_error_list,
                 color='green',
                 label="Visible",
                 linewidth=2,
                 linestyle="None")

    # using math_eq  to seperate data by galaxy to make an average
    x_v_data_seperated = math_eq.split_data(xV_data_list, length)
    x_ir_data_seperated = math_eq.split_data(xIR_data_list, length)
    x_ir_data_err_seperated = math_eq.split_data(xIR_data_error_list, length)
    x_v_data_err_seperated = math_eq.split_data(xV_data_error_list, length)

    average_distanceV = []
    average_distanceIR = []
    average_distanceV_err = []
    average_distanceIR_err = []
    seperated = []
    avg_distance_err = []
    
    # for loop that adds the average distance all into infrared and visible data
    # lists and their errors of the seperated data into their respective galaxies
    for i in range(len(x_v_data_seperated)):
        average_distanceV.append(math_eq.inverse_mean(x_v_data_seperated[i],
                                                      x_v_data_err_seperated[i]))
        average_distanceV_err.append(math_eq.inverse_mean_error(x_v_data_err_seperated[i]))
        average_distanceIR.append(math_eq.inverse_mean(x_ir_data_seperated[i],
                                                       x_ir_data_err_seperated[i]))
        average_distanceIR_err.append(math_eq.inverse_mean_error(x_ir_data_err_seperated[i]))
        seperated.append(np.mean(x_v_data_seperated[i]))

    avg_distance = []

    # for loop that addes the infra red and visible data to one list and same
    # with error so that they can be plotted as a single point that is averaged in
    # their respective galaxies
    for i in range(len(average_distanceIR)):
        avg_distance.append(math_eq.inverse_mean([average_distanceIR[i],
                                                  average_distanceV[i]],
                                                 [average_distanceIR_err[i],
                                                  average_distanceV_err[i]]))
        avg_distance_err.append(math_eq.inverse_mean_error([average_distanceIR_err[i],
                                                            average_distanceV_err[i]]))

    plt.legend()
    plt.show()

    # plotting error bar
    plt.errorbar(avg_distance,
                 velocity,
                 xerr=avg_distance_err,
                 linestyle="None",
                 marker="x")

    # finding the fit and covariance of hubble's law
    fit, covar = opt.curve_fit(math_eq.hubbleLaw,
                               avg_distance,
                               velocity,
                               sigma=avg_distance_err)

    # return row of 100 evenly spaces data of the minimum and max of average distance
    x = np.linspace(min(avg_distance),
                    max(avg_distance),
                    100)
    # plotting fit
    plt.plot(x,
             x * fit,
             color='orange',
             label="Hubble's constant")

    # error= square root of variance
    error = np.sqrt(covar)

    # plotting labels
    plt.xlabel("Average distance per galaxy Mpc")
    plt.ylabel("Recessional Velocity km/s")
    plt.title("Hubble's Constant")

    print("Hubble's constant is", fit, "±", error, " km/s/Mpc")
    # calcualting age of universe in billions of years
    age = (1 / (fit / 3.24078 * 10 ** -20)) / (3.154 * 10 ** 7) / (10 ** 10)

    # proportional error - 1/a doesnt change the proportionality of the error
    age_err = (error / fit) * age

    print("The age of the universe is ", age, "±", age_err, " billion years")

    # calculating chi using math_eq
    chi = math_eq.chisq(avg_distance, np.array(velocity) / fit,
                        np.array(avg_distance_err))

    print("chisq for Hubble's Constant is ", round(chi, 3))
    # calculating reduced chi
    print("Reduced infrared Chi Square is ", round(chi / (len(avg_distance) - 1), 3))
    plt.legend()
    plt.show()

# running function
if __name__ == '__main__':
    get_Hubbles_constant()