# Cox-Ingersoll-Ross Project

### Calibrating and Projecting the Cox-Ingersoll-Ross Process

Before running this script, you must specify your 'working directory'. This is the location where any output will be stored.

*For Mac: 'MISC'->'Change Working Directory...' (And then specify location)*

*For PC: Exectue the following line of code: setwd('location'). Excute by "CNTRL+R".*

The only thing you must specify (besides working directory) is the 'ncap' - or number of periods you wish to forecast - and the 'sim' - the number of simulations you wish to run

#### Running the file...

*For Mac: Highlight all of the lines of code and then 'Edit'->'Execute'*

*For PC: 'Edit'->'Run All'*

**Note**: You will then be prompted to choose your data file in the browser. Once you specify the file, R will then continue to execute the remaining lines of code. Your data file must be saved as a ".csv" and interest rates you wish to calibrate/forecast must be named 'Data' in the column header. Also, you must have a "Counter" column from 1 to as many data entries you have.

Your ouput will be saved as a ".png" and ".txt" file in your working directory.


### Results

Results from CIR_Calibration_Projection.R.

![](https://raw.githubusercontent.com/MGallow/CIR/master/CIR_Calibration_Projections.png)

Results from CIR_Pricing_Projections(monthly).R.

![](https://raw.githubusercontent.com/MGallow/CIR/master/CIR_Pricing_Projections(monthly).png)
