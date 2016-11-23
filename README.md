# icesProxy
ICES methods for setting proxy MSY reference points for stocks in categories 3 and 4

ICES provides advice on fishing opportunities for more than 200 stocks and strives to base this advice on stock status relative to Maximum Sustainable Yield (MSY) objectives. More than 30% of stocks ICES gives advice on have insufficient data and knowledge to conduct a full analytical assessment of their state and exploitation ([ICES category 3 and 4 stocks](http://www.ices.dk/sites/pub/Publication%20Reports/Expert%20Group%20Report/acom/2012/ADHOC/DLS%20Guidance%20Report%202012.pdf)). Methods to provide reference points for such stocks are nascent and are new to most in the scientific community. ICES began to develop methods in 2015 ([WKLIFE V](http://ices.dk/sites/pub/Publication%20Reports/Expert%20Group%20Report/acom/2015/WKLIFEV/wklifeV_2015.pdf)) and applied these methods to a small group of stocks via a workshop ([WKProxy](http://www.ices.dk/sites/pub/Publication%20Reports/Expert%20Group%20Report/acom/2015/WKProxy/01%20WKProxy%20Report.pdf)) to classify the status of category 3 and 4 stocks relative to MSY reference points for 2016. In 2017, ICES will apply these methods for all category 3 and 4 stocks with assessments that should be updated in 2017.

## Organization
ICES uses two major categories of methods for setting proxy MSY reference points: (1) length-based methods, and (2) catch and survey-based methods.

The ICES Technical Guidance on Proxy Reference Point Methodologies will provide in-depth information on the methods and how to interpret the results. ADD LINK

The DTU instructors developed www.stockassessment.org, which houses SAM and SPiCT 
Participants will be given access to this web interface for running SPiCT and accessing examples. The SPiCT source code is open source and distributed as an R package, so it is also possible to run SPiCT with or without using the web interface.


*	Overview
 *	ICES approach to proxy reference points
 *	Orientation to the ICES Technical Guidance
 *	Methods introduction
 *	ICES Advice Rule for proxy reference points (how these values will be used in the advice)

###	Length-based methods

*	Mean-length Z (Gedamke and Hoenig (2006))
 *	Data needed = time series of mean length (>full selection) in catch, length of full selection (estimated from catch length frequency data), growth parameters
 *	Assumptions = recruitment and fishery selection constant over time, fishery selection is steep and asymptotic
 *	Output = total mortality Z
 *	Data format required
 *	Running the model
 *	Interpreting the results

*	Length-based spawner per recruit (LB-SPR)
 *	Data needed = proportion-at-length in catch (in a single year) and biological parameters: L<sub>inf</sub> and M/K
 *	Assumptions = recruitment, fishery selection and F constant over time, fishery selection asymptotic
 *	Output = estimates F/M and compares with reference points from Per-Recruit analysis (requires weight and maturity parameters). 
 *	Data format required
 *	Running the model
 *	Interpreting the results

*	Length-based indicators (LBI)
 *	Data needed = proportion-at-length in catch (single year), length at first capture in fishery (estimated from length frequency data), biological parameters (L<sub>inf</sub> and length-at-maturity)
 *	Assumptions = recruitment, fishery selection and F constant over time, fishery selection asymptotic
 *	Output = Range of indicators and their expected values when exploitation is consistent with conservation of large and immature individuals, optimal yield and MSY
 *	Data format required
 *	Running the model
 *	Interpreting the results

*	S6
 *	Data needed
 *	Assumptions
 *	Provides 
 *	Data format required
 *	Running the model
 *	Interpreting the results

###	Biomass dynamics models

*	CMSY
 *  Data needed = time series of catch
 *  Assumptions = many
 *  Provides F/F<sub>MSY</sub> and B/MSY B<sub>trigger</sub>
 *  Data format required
 *  Running the model
 *	Interpreting the results

*	SPiCT
 *	Data needed = time series of catch (by weight) and biomass index
 *	Assumptions = constant fishery selectivity over time, biomass index proportional to exploitable stock biomass
 *	Provides F/F<sub>MSY</sub> and B/MSY B<sub>trigger</sub>
 *	Data format required
 *	Running the model
 *	Interpreting the results
