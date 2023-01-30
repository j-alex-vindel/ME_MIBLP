import sys
import os
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'MY_MODULES'))) 

from Ob_Met_Net import Metabolic_Network
import pandas as pd

yeast = 'Yeast'
folder_yeast = "../B_Yeast/Data/Yeast"


# =========== Access Data ====================

S = pd.read_csv(f'{folder_yeast}/yeast_S.csv',header=None).values
rxn = pd.read_csv(f'{folder_yeast}/yeast_rxns.csv',header=None)[0].tolist()
met = pd.read_csv(f'{folder_yeast}/yeast_mets.csv',header=None)[0].tolist()
rnames = pd.read_csv(f'{folder_yeast}/yeast_rxnNames.csv',header=None)[0].tolist()
LB = pd.read_csv(f'{folder_yeast}/LB_biolow_y.txt',header=None)[0].tolist()
UB = pd.read_csv(f'{folder_yeast}/UB_y.txt',header=None)[0].tolist()

biomas = rxn.index('r_2110')
growth = rnames.index('growth')
chemical = rxn.index('r_1761')


essential = [1542,1545,1546,1547,1548,1549,1550,1551,1552,1553,1554,1563,1564,1565,1566,
            1571,1572,1577,1580,1581,1586,1589,1598,1604,1621,1625,1627,1629,1630,1631,
            1634,1639,1641,1643,1648,1649,1650,1651,1654,1663,1671,1672,1683,1687,1702,
            1705,1706,1709,1710,1711,1712,1714,1715,1716,1718,1727,1730,1749,1753,1757,
            1761,1764,1765,1767,1788,1791,1792,1793,1798,1799,1800,1806,1807,1808,1810,
            1814,1815,1818,1820,1831,1833,1840,1842,1846,1860,1861,1864,1865,1866,1869,
            1872,1874,1877,1878,1879,1880,1882,1885,1888,1890,1892,1895,1896,1898,1899,
            1901,1902,1903,1905,1908,1910,1911,1912,1913,1914,1915,1916,1917,1930,1932,
            1946,1948,1950,1951,1961,1966,1967,1983,1984,1986,1988,1991,1992,1993,1998,
            1999,2000,2004,2008,2018,2019,2023,2027,2032,2037,2042,2043,2045,2048,2050,
            2051,2054,2055,2057,2059,2060,2061,2065,2066,2067,2072,2075,2082,2089,2090,
            2091,2099,2101,2103,2105,2109,2110]

T = [i-1 for i in essential]

knockout = [i for i in range(len(rxn)) if i not in T]

# =========== Metabolic Network Object ====================

MetNet_Yeast = Metabolic_Network(S=S,LB=LB,UB=UB,Rxn=rxn,Met=met,Name=yeast,KO=knockout,biomass= biomas, chemical=chemical)

