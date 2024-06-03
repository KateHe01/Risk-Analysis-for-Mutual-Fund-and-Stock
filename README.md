# Risk Analysis for Mutual Fund and Stock

This project performs a risk analysis on a mutual fund (Fund X) and a stock (Stock G) using R. The analysis involves estimating Value at Risk (VaR) and Expected Shortfall (ES) for both investments.

## Project Structure

The project is divided into three main parts:
1. Analysis of Fund X (Question 1, 2, and 3)
2. Analysis of Stock G (Question 4)
3. Rolling VaR Analysis (Question 5)

### Part 1: Analysis of Fund X

#### Data
- `logret_X.xts`: Daily log returns of Fund X from 2013-01-31 to 2018-01-30.
- `holdings_X`: Holdings of Fund X as of 2017-10-31.

#### Steps

1. **Loading Data**: Load the data objects into the R workspace.
2. **Analyzing Log Returns**:
    - Calculate basic statistics (mean, variance) of the log returns.
    - Determine the most appropriate statistical model for the returns.
    - Estimate the 95% VaR and ES for the next five trading days.
3. **Estimating VaR and ES**: Use historical simulation to estimate VaR and ES.
4. **Implied Volatility Calculation**:
    - Calculate the implied volatilities of options using Black's pricing model.
    - Check for the presence of volatility skew.
5. **Scenario Analysis**: 
    - Analyze the impact of increased implied volatility on the fund’s value.
    - Determine how much the implied volatilities need to rise to reduce the fund value by 50%.

### Results for Fund X

- **VaR (95%)**: -0.0469
- **ES (95%)**: -0.0851
- **Scenario Analysis**:
  - **Increase by 50%**: Total Fund Value = 544,193,636
  - **Increase by 100%**: Total Fund Value = 401,948,788
  - **Increase by 200%**: Total Fund Value = 6,104,915
  - **Increase by 300%**: Total Fund Value = -460,537,605
- **Implied Volatility Increase Needed to Reduce Fund Value by 50%**: 125%

### Part 2: Analysis of Stock G

#### Data
- `logret_G.xts`: Daily log returns of Stock G from 2002-03-04 to 2021-06-30.

#### Steps

1. **Loading Data**: Load the data objects into the R workspace.
2. **Estimating VaR using GARCH Model**:
    - Fit a GARCH-t model to the return series.
    - Use the rolling forecast to estimate 1-day ahead VaR at the 95% confidence level.
3. **Calculating Margin Requirements**:
    - Calculate the margin requirement for a $1,000,000 position in Stock G using the estimated VaR.

### Results for Stock G

- **Margin Requirements**: Calculated based on the VaR estimates for a $1,000,000 position.

### Part 3: Rolling VaR Analysis

#### Steps

1. **Rolling VaR Forecast**:
    - Perform a rolling forecast using a GARCH-t model to estimate 1-day ahead VaR.
2. **Plotting the Results**:
    - Plot the 1-day 95% VaR alongside the actual returns to visualize breaches.
3. **Counting VaR Breaches**:
    - Count the number of days where actual returns fall below the 95% VaR from the rolling forecast.

### Results for Rolling VaR Analysis

- **Plot of 1-day 95% VaR from Dec 31, 2020 to Jun 30, 2021**:
- **Number of VaR breaches**: 7 out of 124 days

### Conclusion

This project demonstrates the application of risk management techniques to assess the risk of a mutual fund and a stock using historical data and statistical models. The analysis includes statistical modeling, scenario analysis, and rolling VaR forecasts to provide a comprehensive view of the risks associated with the investments.
