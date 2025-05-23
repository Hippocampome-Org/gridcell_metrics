# Gridcell Metrics Reporting

Current version: 1.0.3

## Related articles
<br>Some additional details and analyses run with this software can be found in this article.
<br>\[1\] R. G., R., Ascoli, G. A., Sutton, N. M., & Dannenberg, H. (2024). Spatial periodicity in grid cell firing is explained by a neural sequence code of 2-D trajectories. eLife, 13. [https://doi.org/10.7554/eLife.96627.1](https://doi.org/10.7554/eLife.96627.1)

<br>An article describing methods in and usage of the software is available here:
<br>\[2\] Sutton, N. M., Gutiérrez-Guzmán, B. E., Dannenberg, H., & Ascoli, G. A. (2025). Automated Measurement of Grid Cell Firing Characteristics. Algorithms, 18(3), 139. [https://doi.org/10.3390/a18030139](https://doi.org/10.3390/a18030139)

## Requirements
<br>If users are launching this software from MATLAB then the DBSCAN clustering algorithm from MATLAB’s (mathworks.com) Statistics and Machine Learning Toolbox is needed. A user should install that toolbox in that case. If a user uses the runtime version of the software, which can be downloaded under the releases section, then the toolbox is not needed to be installed.

## Example Usage Instructions
<br>See [instructions](https://hco-dev-docs.readthedocs.io/en/latest/gridcell_metrics/usage_instruct.html) for further details

## General Usage
<br>src/run_gridcell_metrics.m should be run in MATLAB for running the software from source code. run_gridcell_metrics_with_config.bat (Windows) or run_gridcell_metrics_with_config.sh (Linux) should be used to run the runtime versions of this software. The runtime versions of this software are available in the releases section of Github.

## Recreating Results Reported in (R. G. et al., 2024)
<br>How to recreate the results that found 31% is a top-performing non-field value filter threshold is described in [recreating results](https://hco-dev-docs.readthedocs.io/en/latest/gridcell_metrics/recreating_results.html).

## Additional Info.
<br>See the [documentation](https://hco-dev-docs.readthedocs.io/en/latest/gridcell_metrics/overview.html) site for further info.