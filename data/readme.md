# SolarChem Data readme


## Astromat

`astromat/` contains files from the [Astromaterials Data System](https://www.astromat.org)


| file | description |
| :--- | :---------- |
| `astromat-6_27_2024-eventId-2294.csv` | A csv of exported AstroMat data, accessed June 27, 2024. |
| `fastromat.jl` | Julia script to update `fastromat.jls`|
| `fastromat.jls` | A (Julia) serialized file of API-accessed Astromat data. Used with the `SolarChem.fastromat` function.
| `taxon-num-name.txt` | Taxon names corresponding to taxon-num IDs from Astromat API. |

## Solar Twins

`Bedell2018-solar-twins.csv` contains compiled data from Bedell + 2018 (*ApJ*, doi:[10.3847/1538-4357/aad908](https://doi.org/10.3847/1538-4357/aad908)). 