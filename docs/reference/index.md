# Package index

## Data setup

Functions to set up the Seurat object for network analysis

- [`SetupForWGCNA()`](https://smorabit.github.io/hdWGCNA/reference/SetupForWGCNA.md)
  : SetupForWGCNA
- [`SelectNetworkGenes()`](https://smorabit.github.io/hdWGCNA/reference/SelectNetworkGenes.md)
  : SelectNetworkGenes
- [`FindMajorIsoforms()`](https://smorabit.github.io/hdWGCNA/reference/FindMajorIsoforms.md)
  : FindMajorIsoforms

## Metacells and metaspots

Functions for constructing metacells from single-cell data and metaspots
from ST data

- [`MetacellsByGroups()`](https://smorabit.github.io/hdWGCNA/reference/MetacellsByGroups.md)
  : MetacellsByGroups
- [`MetaspotsByGroups()`](https://smorabit.github.io/hdWGCNA/reference/MetaspotsByGroups.md)
  : MetaspotsByGroups

## Network Analysis

Functions for constructing the co-expression network

- [`TestSoftPowers()`](https://smorabit.github.io/hdWGCNA/reference/TestSoftPowers.md)
  : TestSoftPowers
- [`TestSoftPowersConsensus()`](https://smorabit.github.io/hdWGCNA/reference/TestSoftPowersConsensus.md)
  : TestSoftPowersConsensus
- [`PlotSoftPowers()`](https://smorabit.github.io/hdWGCNA/reference/PlotSoftPowers.md)
  : PlotSoftPowers
- [`ConstructNetwork()`](https://smorabit.github.io/hdWGCNA/reference/ConstructNetwork.md)
  : ConstructNetwork
- [`ModuleEigengenes()`](https://smorabit.github.io/hdWGCNA/reference/ModuleEigengenes.md)
  : ModuleEigengenes
- [`ModuleConnectivity()`](https://smorabit.github.io/hdWGCNA/reference/ModuleConnectivity.md)
  : ModuleConnectivity

## Network Visualization

Functions for visualizing the co-expression network

- [`PlotDendrogram()`](https://smorabit.github.io/hdWGCNA/reference/PlotDendrogram.md)
  : PlotDendrogram
- [`ModuleNetworkPlot()`](https://smorabit.github.io/hdWGCNA/reference/ModuleNetworkPlot.md)
  : ModuleNetworkPlot
- [`HubGeneNetworkPlot()`](https://smorabit.github.io/hdWGCNA/reference/HubGeneNetworkPlot.md)
  : HubGeneNetworkPlot
- [`RunModuleUMAP()`](https://smorabit.github.io/hdWGCNA/reference/RunModuleUMAP.md)
  : RunModuleUMAP
- [`ModuleUMAPPlot()`](https://smorabit.github.io/hdWGCNA/reference/ModuleUMAPPlot.md)
  : ModuleUMAPPlot

## Differential analysis

Functions for differential module eigengene analysis

- [`FindDMEs()`](https://smorabit.github.io/hdWGCNA/reference/FindDMEs.md)
  : FindDMEs
- [`FindAllDMEs()`](https://smorabit.github.io/hdWGCNA/reference/FindAllDMEs.md)
  : FindAllDMEs
- [`PlotDMEsVolcano()`](https://smorabit.github.io/hdWGCNA/reference/PlotDMEsVolcano.md)
  : PlotDMEsVolcano
- [`PlotDMEsLollipop()`](https://smorabit.github.io/hdWGCNA/reference/PlotDMEsLollipop.md)
  : PlotDMEsLollipop

## Enrichment Analysis

Functions for Enrichr analysis

- [`RunEnrichr()`](https://smorabit.github.io/hdWGCNA/reference/RunEnrichr.md)
  : RunEnrichr
- [`EnrichrBarPlot()`](https://smorabit.github.io/hdWGCNA/reference/EnrichrBarPlot.md)
  : EnrichrBarPlot
- [`EnrichrDotPlot()`](https://smorabit.github.io/hdWGCNA/reference/EnrichrDotPlot.md)
  : EnrichrDotPlot
- [`OverlapModulesDEGs()`](https://smorabit.github.io/hdWGCNA/reference/OverlapModulesDEGs.md)
  : OverlapModulesDEGs
- [`OverlapBarPlot()`](https://smorabit.github.io/hdWGCNA/reference/OverlapBarPlot.md)
  : OverlapBarPlot
- [`OverlapDotPlot()`](https://smorabit.github.io/hdWGCNA/reference/OverlapDotPlot.md)
  : OverlapDotPlot

## Plotting

Functions for generating plots with hdWGCNA

- [`ModuleFeaturePlot()`](https://smorabit.github.io/hdWGCNA/reference/ModuleFeaturePlot.md)
  : ModuleFeaturePlot
- [`ModuleRadarPlot()`](https://smorabit.github.io/hdWGCNA/reference/ModuleRadarPlot.md)
  : ModuleRadarPlot
- [`PlotModuleTrajectory()`](https://smorabit.github.io/hdWGCNA/reference/PlotModuleTrajectory.md)
  : PlotModuleTrajectory

## Transcription Factor Regulatory Networks

- [`MotifScan()`](https://smorabit.github.io/hdWGCNA/reference/MotifScan.md)
  : MotifScan
- [`ConstructTFNetwork()`](https://smorabit.github.io/hdWGCNA/reference/ConstructTFNetwork.md)
  : ConstructTFNetwork
- [`AssignTFRegulons()`](https://smorabit.github.io/hdWGCNA/reference/AssignTFRegulons.md)
  : AssignTFRegulons
- [`GetTFTargetGenes()`](https://smorabit.github.io/hdWGCNA/reference/GetTFTargetGenes.md)
  : GetTFTargetGenes
- [`TFNetworkPlot()`](https://smorabit.github.io/hdWGCNA/reference/TFNetworkPlot.md)
  : TFNetworkPlot
- [`RegulonScores()`](https://smorabit.github.io/hdWGCNA/reference/RegulonScores.md)
  : RegulonScores
- [`RegulonBarPlot()`](https://smorabit.github.io/hdWGCNA/reference/RegulonBarPlot.md)
  : RegulonBarPlot
- [`RunEnrichrRegulons()`](https://smorabit.github.io/hdWGCNA/reference/RunEnrichrRegulons.md)
  : RunEnrichrRegulons
- [`ModuleRegulatoryNetwork()`](https://smorabit.github.io/hdWGCNA/reference/ModuleRegulatoryNetwork.md)
  : ModuleRegulatoryNetwork
- [`ModuleRegulatoryNetworkPlot()`](https://smorabit.github.io/hdWGCNA/reference/ModuleRegulatoryNetworkPlot.md)
  : ModuleRegulatoryNetworkPlot
- [`ModuleRegulatoryHeatmap()`](https://smorabit.github.io/hdWGCNA/reference/ModuleRegulatoryHeatmap.md)
  : ModuleRegulatoryHeatmap
- [`FindDifferentialRegulons()`](https://smorabit.github.io/hdWGCNA/reference/FindDifferentialRegulons.md)
  : FindDifferentialRegulons
- [`PlotDifferentialRegulons()`](https://smorabit.github.io/hdWGCNA/reference/PlotDifferentialRegulons.md)
  : PlotDifferentialRegulons

## Module Preservation

Functions for performing module preservation analysis

- [`ProjectModules()`](https://smorabit.github.io/hdWGCNA/reference/ProjectModules.md)
  : ProjectModules
- [`ModulePreservation()`](https://smorabit.github.io/hdWGCNA/reference/ModulePreservation.md)
  : ModulePreservation
- [`ModulePreservationNetRep()`](https://smorabit.github.io/hdWGCNA/reference/ModulePreservationNetRep.md)
  : ModulePreservationNetRep
- [`PlotModulePreservation()`](https://smorabit.github.io/hdWGCNA/reference/PlotModulePreservation.md)
  : PlotModulePreservation
- [`PlotModulePreservationLollipop()`](https://smorabit.github.io/hdWGCNA/reference/PlotModulePreservationLollipop.md)
  : PlotModulePreservationLollipop
- [`ModuleTopologyHeatmap()`](https://smorabit.github.io/hdWGCNA/reference/ModuleTopologyHeatmap.md)
  : ModuleTopologyHeatmap
- [`ModuleTopologyBarplot()`](https://smorabit.github.io/hdWGCNA/reference/ModuleTopologyBarplot.md)
  : ModuleTopologyBarplot

## Module Trait Correlation

Functions for performing module trait correlation analysis

- [`ModuleTraitCorrelation()`](https://smorabit.github.io/hdWGCNA/reference/ModuleTraitCorrelation.md)
  : Module-Trait Correlation
- [`PlotModuleTraitCorrelation()`](https://smorabit.github.io/hdWGCNA/reference/PlotModuleTraitCorrelation.md)
  : PlotModuleTraitCorrelation

## Transcription Factor Networks

- [`MotifScan()`](https://smorabit.github.io/hdWGCNA/reference/MotifScan.md)
  : MotifScan

## Getters and setters

Functions to retrieve and set the values for various attributes

- [`SetActiveWGCNA()`](https://smorabit.github.io/hdWGCNA/reference/SetActiveWGCNA.md)
  : SetActiveWGCNA

- [`GetActiveWGCNA()`](https://smorabit.github.io/hdWGCNA/reference/GetActiveWGCNA.md)
  : GetActiveWGCNA

- [`SetMetacellObject()`](https://smorabit.github.io/hdWGCNA/reference/SetMetacellObject.md)
  : SetMetacellObject

- [`GetMetacellObject()`](https://smorabit.github.io/hdWGCNA/reference/GetMetacellObject.md)
  : GetMetacellObject

- [`SetWGCNAGenes()`](https://smorabit.github.io/hdWGCNA/reference/SetWGCNAGenes.md)
  : SetWGCNAGenes

- [`GetWGCNAGenes()`](https://smorabit.github.io/hdWGCNA/reference/GetWGCNAGenes.md)
  : GetWGCNAGenes

- [`SetDatExpr()`](https://smorabit.github.io/hdWGCNA/reference/SetDatExpr.md)
  : Set Expression Data (Standard WGCNA)

- [`GetDatExpr()`](https://smorabit.github.io/hdWGCNA/reference/GetDatExpr.md)
  : GetDatExpr

- [`SetMultiExpr()`](https://smorabit.github.io/hdWGCNA/reference/SetMultiExpr.md)
  :

  This function prepares the expression data for Consensus WGCNA
  analysis. It populates the `multiExpr` slot in the Seurat object,
  which contains a list of expression matrices (one for each consensus
  group, e.g., dataset, sample, condition).

- [`GetMultiExpr()`](https://smorabit.github.io/hdWGCNA/reference/GetMultiExpr.md)
  : GetMultiExpr

- [`SetWGCNAParams()`](https://smorabit.github.io/hdWGCNA/reference/SetWGCNAParams.md)
  : SetWGCNAParams

- [`GetWGCNAParams()`](https://smorabit.github.io/hdWGCNA/reference/GetWGCNAParams.md)
  : GetWGCNAParams

- [`SetPowerTable()`](https://smorabit.github.io/hdWGCNA/reference/SetPowerTable.md)
  : SetPowerTable

- [`GetPowerTable()`](https://smorabit.github.io/hdWGCNA/reference/GetPowerTable.md)
  : GetPowerTable

- [`SetNetworkData()`](https://smorabit.github.io/hdWGCNA/reference/SetNetworkData.md)
  : SetNetworkData

- [`GetNetworkData()`](https://smorabit.github.io/hdWGCNA/reference/GetNetworkData.md)
  : GetNetworkData

- [`SetModules()`](https://smorabit.github.io/hdWGCNA/reference/SetModules.md)
  : SetModules

- [`GetModules()`](https://smorabit.github.io/hdWGCNA/reference/GetModules.md)
  : GetModules

- [`GetHubGenes()`](https://smorabit.github.io/hdWGCNA/reference/GetHubGenes.md)
  : GetHubGenes

- [`SetMEs()`](https://smorabit.github.io/hdWGCNA/reference/SetMEs.md) :
  SetMEs

- [`GetMEs()`](https://smorabit.github.io/hdWGCNA/reference/GetMEs.md) :
  GetMEs

- [`SetMELoadings()`](https://smorabit.github.io/hdWGCNA/reference/SetMELoadings.md)
  : SetMELoadings

- [`GetMELoadings()`](https://smorabit.github.io/hdWGCNA/reference/GetMELoadings.md)
  : GetMELoadings

- [`SetEnrichrTable()`](https://smorabit.github.io/hdWGCNA/reference/SetEnrichrTable.md)
  : SetEnrichrTable

- [`GetEnrichrTable()`](https://smorabit.github.io/hdWGCNA/reference/GetEnrichrTable.md)
  : GetEnrichrTable

- [`SetModuleScores()`](https://smorabit.github.io/hdWGCNA/reference/SetModuleScores.md)
  : SetModuleScores

- [`GetModuleScores()`](https://smorabit.github.io/hdWGCNA/reference/GetModuleScores.md)
  : GetModuleScores

- [`GetTOM()`](https://smorabit.github.io/hdWGCNA/reference/GetTOM.md) :
  GetTOM

- [`SetModuleUMAP()`](https://smorabit.github.io/hdWGCNA/reference/SetModuleUMAP.md)
  : SetModuleUMAP

- [`GetModuleUMAP()`](https://smorabit.github.io/hdWGCNA/reference/GetModuleUMAP.md)
  : GetModuleUMAP

- [`SetModuleTraitCorrelation()`](https://smorabit.github.io/hdWGCNA/reference/SetModuleTraitCorrelation.md)
  : SetModuleTraitCorrelation

- [`GetModuleTraitCorrelation()`](https://smorabit.github.io/hdWGCNA/reference/GetModuleTraitCorrelation.md)
  : GetModuleTraitCorrelation

- [`SetModulePreservation()`](https://smorabit.github.io/hdWGCNA/reference/SetModulePreservation.md)
  : SetModulePreservation

- [`GetModulePreservation()`](https://smorabit.github.io/hdWGCNA/reference/GetModulePreservation.md)
  : GetModulePreservation

- [`SetMotifs()`](https://smorabit.github.io/hdWGCNA/reference/SetMotifs.md)
  : SetMotifs

- [`GetMotifs()`](https://smorabit.github.io/hdWGCNA/reference/GetMotifs.md)
  : GetMotifs

- [`SetMotifMatrix()`](https://smorabit.github.io/hdWGCNA/reference/SetMotifMatrix.md)
  : SetMotifMatrix

- [`GetMotifMatrix()`](https://smorabit.github.io/hdWGCNA/reference/GetMotifMatrix.md)
  : GetMotifMatrix

- [`SetPFMList()`](https://smorabit.github.io/hdWGCNA/reference/SetPFMList.md)
  : SetPFMList

- [`GetPFMList()`](https://smorabit.github.io/hdWGCNA/reference/GetPFMList.md)
  : GetPFMList

- [`SetMotifTargets()`](https://smorabit.github.io/hdWGCNA/reference/SetMotifTargets.md)
  : SetMotifTargets

- [`GetMotifTargets()`](https://smorabit.github.io/hdWGCNA/reference/GetMotifTargets.md)
  : GetMotifTargets

- [`SetMotifOverlap()`](https://smorabit.github.io/hdWGCNA/reference/SetMotifOverlap.md)
  : SetMotifOverlap

- [`GetMotifOverlap()`](https://smorabit.github.io/hdWGCNA/reference/GetMotifOverlap.md)
  : GetMotifOverlap

- [`SetMotifScores()`](https://smorabit.github.io/hdWGCNA/reference/SetMotifScores.md)
  : SetMotifScores

- [`GetMotifScores()`](https://smorabit.github.io/hdWGCNA/reference/GetMotifScores.md)
  : GetMotifScores

- [`SetDegrees()`](https://smorabit.github.io/hdWGCNA/reference/SetDegrees.md)
  : SetDegrees

- [`GetDegrees()`](https://smorabit.github.io/hdWGCNA/reference/GetDegrees.md)
  : GetDegrees

- [`GetTFRegulons()`](https://smorabit.github.io/hdWGCNA/reference/GetTFRegulons.md)
  : GetTFRegulons

- [`SetTFRegulons()`](https://smorabit.github.io/hdWGCNA/reference/SetTFRegulons.md)
  : SetTFRegulons

- [`GetTFEval()`](https://smorabit.github.io/hdWGCNA/reference/GetTFEval.md)
  : GetTFEval

- [`SetTFEval()`](https://smorabit.github.io/hdWGCNA/reference/SetTFEval.md)
  : SetTFEval

- [`GetRegulonScores()`](https://smorabit.github.io/hdWGCNA/reference/GetRegulonScores.md)
  : GetRegulonScores

- [`SetRegulonScores()`](https://smorabit.github.io/hdWGCNA/reference/SetRegulonScores.md)
  : SetRegulonScores

- [`GetEnrichrRegulonTable()`](https://smorabit.github.io/hdWGCNA/reference/GetEnrichrRegulonTable.md)
  : GetEnrichrRegulonTable

- [`SetEnrichrRegulonTable()`](https://smorabit.github.io/hdWGCNA/reference/SetEnrichrRegulonTable.md)
  : SetEnrichRegulonTable

- [`GetMetacellParams()`](https://smorabit.github.io/hdWGCNA/reference/GetMetacellParams.md)
  : GetMetacellParams

- [`SetMetacellParams()`](https://smorabit.github.io/hdWGCNA/reference/SetMetacellParams.md)
  : SetMetacellParams

- [`SetTFNetwork()`](https://smorabit.github.io/hdWGCNA/reference/SetTFNetwork.md)
  : SetTFNetwork

- [`GetTFNetwork()`](https://smorabit.github.io/hdWGCNA/reference/GetTFNetwork.md)
  : GetTFNetwork

## Seurat wrappers

Wrapper functions to run Seurat commands on the metacell data

- [`NormalizeMetacells()`](https://smorabit.github.io/hdWGCNA/reference/NormalizeMetacells.md)
  : NormalizeMetacells
- [`ScaleMetacells()`](https://smorabit.github.io/hdWGCNA/reference/ScaleMetacells.md)
  : ScaleMetacells
- [`RunPCAMetacells()`](https://smorabit.github.io/hdWGCNA/reference/RunPCAMetacells.md)
  : RunPCAMetacells
- [`RunHarmonyMetacells()`](https://smorabit.github.io/hdWGCNA/reference/RunHarmonyMetacells.md)
  : RunHarmonyMetacells
- [`RunUMAPMetacells()`](https://smorabit.github.io/hdWGCNA/reference/RunUMAPMetacells.md)
  : RunUMAPMetacells
- [`DimPlotMetacells()`](https://smorabit.github.io/hdWGCNA/reference/DimPlotMetacells.md)
  : DimPlotMetacells

## Customization

Customization functions

- [`ResetModuleColors()`](https://smorabit.github.io/hdWGCNA/reference/ResetModuleColors.md)
  : ResetModuleColors
- [`ResetModuleNames()`](https://smorabit.github.io/hdWGCNA/reference/ResetModuleNames.md)
  : ResetModuleNames

## Other

Other functions

- [`AggregatePseudobulk()`](https://smorabit.github.io/hdWGCNA/reference/AggregatePseudobulk.md)
  : AggregatePseudobulk
- [`NormalizeCounts()`](https://smorabit.github.io/hdWGCNA/reference/NormalizeCounts.md)
  : NormalizeCounts
- [`ReassignModules()`](https://smorabit.github.io/hdWGCNA/reference/ReassignModules.md)
  : ReassignModules
- [`GetActiveWGCNAName()`](https://smorabit.github.io/hdWGCNA/reference/GetActiveWGCNAName.md)
  : GetActiveWGCNAName
