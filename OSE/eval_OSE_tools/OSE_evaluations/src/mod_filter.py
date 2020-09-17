import xarray as xr
import numpy as np
import scipy.signal
from scipy.signal import butter, filtfilt
from functools import partial


from bokeh.layouts import gridplot
from bokeh.plotting import figure, show, output_file
from bokeh.models import HoverTool, LinearAxis, Range1d
from bokeh.models.formatters import DatetimeTickFormatter


def apply_lowpass_filter(data_raw, time, cutoff_freq):
    """

    :param data_raw:
    :param time:
    :param cutoff_freq:
    :return:
    """

    display = False
    # Buterworth filter
    order = 5  # Filter order
    b, a = scipy.signal.butter(order, cutoff_freq[0], output='ba', btype='lowpass')

    # data size must be greater than filter size
    if data_raw.size >  3*b.size:
        
        arr = xr.DataArray(data_raw, coords=[time], dims=['time'])

        # Apply the filter
        filtered = xr.apply_ufunc(partial(filtfilt, b, a),
                                  arr.chunk(),
                                  dask='parallelized',
                                  output_dtypes=[arr.dtype],
                                  kwargs={'axis': 0}).compute()

        data_filtered = filtered.values

        if display:
            plot_filtering(time, data_raw, data_filtered)
    
        return data_filtered


def apply_highpass_filter(data_raw, time, cutoff_freq):
    """

    :param data_raw:
    :param time:
    :param cutoff_freq:
    :return:
    """

    display = False
    # Buterworth filter
    order = 5  # Filter order
    b, a = scipy.signal.butter(order, cutoff_freq[0], output='ba', btype='highpass')
    
    # data size must be greater than filter size
    if data_raw.size >  3*b.size:

        arr = xr.DataArray(data_raw, coords=[time], dims=['time'])

        # Apply the filter
        filtered = xr.apply_ufunc(partial(filtfilt, b, a),
                                  arr.chunk(),
                                  dask='parallelized',
                                  output_dtypes=[arr.dtype],
                                  kwargs={'axis': 0}).compute()

        data_filtered = filtered.values

        if display:
            plot_filtering(time, data_raw, data_filtered)

        return data_filtered


def apply_bandpass_filter(data_raw, time, cutoff_freq):
    """

    :param data_raw:
    :param time:
    :param cutoff_freq:
    :return:
    """

    display = False
    # Buterworth filter
    order = 5  # Filter order
    b, a = scipy.signal.butter(order, cutoff_freq, output='ba', btype='band')
    
    # data size must be greater than filter size
    if data_raw.size >  3*b.size:
        # print(data_raw.size, b.size)
        arr = xr.DataArray(data_raw, coords=[time], dims=['time'])

        # Apply the filter
        filtered = xr.apply_ufunc(partial(filtfilt, b, a),
                                  arr.chunk(),
                                  dask='parallelized',
                                  output_dtypes=[arr.dtype],
                                  kwargs={'axis': 0}).compute()

        data_filtered = filtered.values
    
        if display:
            plot_filtering(time, data_raw, data_filtered)

        return data_filtered


def plot_filtering(time, raw, filtered):
    """

    :param time:
    :param raw:
    :param filtered:
    :return:
    """
    
    
    hover = HoverTool()
    hover.point_policy='snap_to_data'
    hover.line_policy='nearest'
    TOOLS = ['pan,wheel_zoom,box_zoom']
    
    time = np.arange(0, raw.size)
    p1 = figure(title="Timeseries SSH (raw vs filtered)", tools=TOOLS, plot_width=1000)
    p1.grid.grid_line_alpha=0.3
    p1.xaxis.axis_label = 'Date (Julian days)'
    p1.yaxis.axis_label = 'SSH'
    p1.circle(time,raw, color='lightseagreen')
    p1.line(time, raw, legend='raw', color='lightseagreen')
    p1.circle(time, filtered, color='royalblue')
    p1.line(time, filtered, legend='filtered', color='royalblue')#, line_width=4)

    p1.add_tools(HoverTool(tooltips=[( 'date', '@x'),( 'SSH', '@y m'),],mode='mouse'))

    p1.legend.location = "top_left"
    # output_file("Hss_timeseries.html", title="Hss_timeseries")
    show(p1, plot_width=1000, plot_height=400) 
    
    # Make plots
    # fig = plt.figure()
    # ax1 = fig.add_subplot(211)
    # plt.plot(time, raw, 'b-')
    # plt.plot(time, filtered, 'r-', linewidth=2)
    # plt.ylim(-1,1,0.2)
    # plt.ylabel("SLA (cm)")
    # plt.legend(['Original', 'Filtered'])
    # ax1.axes.get_xaxis().set_visible(False)

    # fig.add_subplot(212)
    # plt.plot(time, raw - filtered, 'b-')
    # plt.ylabel("SLA (cm)")
    # plt.legend(['Residuals'])
    # plt.show()

