# godfa
Go Domain Flow Architecture emulation API server

## Domain Flow Architecture

Domain Flow Architectures are hardware execution models based on the Stillwater distributed data flow engine. 
The Dennis Data Flow Machine (1) is an execution model that theoretically is able to execute along the free schedule of an algorithm.
The free schedule is the inherent parallelism of an algorithm defined by its data dependencies. The model can scale to arbitrary large
computational graphs and arbitrary concurrency.

However, the key disadvantage of the Data Flow Machine architecture is that in physical designs the Content Addressable Memory that 
holds pending instructions, grows with the inherent parallelism of the algorithms, and larger CAMs slow the cycle time of the core 
instruction execution cycle. This has limited data flow machines to just a hundred or so concurrent instructions.

For an important subclass of parallel algorithms, those that exhibit regular data structures and flows, for example, linear algebra,
signal processing, sensor fusion, computer vision, and machine learning algorithms, this regularity can be exploited to create more efficient 
parallel execution by distributing the CAM across the processing elements. This model represents the Stillwater Domain Flow Architecture 
execution model, and in many ways can be seen as the Processor-in-Memory equivalent of the Distributed Memory Machine Architectures that
are currently proposed for exa-scale applications. 

The Domain Flow Architecture is focused on efficient parallel execution to improve computational latency and power consumption by 
recognizing that a distributed machine is bound by the constraints of spacetime. Otherwise stated, information exchange takes time
and energy, and thus it is important that the algorithm designer takes these constraints into account to yield an optimal, energy-efficient
parallel program.


(1) <i>A preliminary architecture for a basic data-flow processor,</i>
Authors:	Jack B. Dennis, David P. Misunas
Project MAC, Massachusetts Institute of Technology
Published in: Proceeding ISCA '75 Proceedings of the 2nd annual symposium on Computer architecture, Pages 126-132 
