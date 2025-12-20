# KmerCO: Tecnica del conteggio dei kmer basata su Bloom Filter

Il progetto nasce dall’esigenza di gestire e analizzare grandi quantità di dati genomici, sfruttando la potenza e l’efficienza dei Bloom Filter e della tecnica KmerCO. Negli anni sono state proposte diverse tecniche, ma KmerCo rappresenta una delle tecniche più all'avanguardia poiché risulta piú precisa ma, allo stesso tempo, richiede una grande quantità di tempo dal punto di vista compitazionale data dalla quantità di confronti che l'algoritmo deve gestire. L’obiettivo è fornire uno strumento quanto più accurato possibile per il conteggio e la classificazione dei k-mer, elementi fondamentali nell’analisi bioinformatica di sequenze di DNA. Verranno proposte due versioni a tal proposito, ovvero la versione Canonica e la versione non Canonica che differiscono tra loro non solo per il costo computazionale impiegato, ma anche per il modo con cui vengono considerate le sequenze genomiche e i loro rispettivi reverse complement. 

## Contesto

Nel campo della bioinformatica, il conteggio dei k-mer (sottosequenze di lunghezza k) è cruciale per molte applicazioni, tra cui l’assemblaggio di genomi, la correzione di errori e la metagenomica. Tuttavia, la gestione di dataset di grandi dimensioni pone sfide significative in termini di memoria e velocità di elaborazione.

## Tecnica Utilizzata: KmerCO

La tecnica KmerCO (K-mer Counting and Organization) si basa sull’utilizzo di Bloom Filter per una rappresentazione compatta e probabilistica dei k-mer. Questo approccio consente di ridurre drasticamente l’occupazione di memoria, mantenendo al contempo un’elevata velocità di accesso e conteggio. Il Bloom Filter permette di verificare rapidamente la presenza di un k-mer, con un tasso di falsi positivi controllato e configurabile.

## Dataset e Test

I test principali sono stati effettuati sul file `balaenoptera.fastq`, scelto per la sua completezza e dimensione significativa (circa 400MB). Questo file rappresenta un caso reale e impegnativo, ideale per valutare le performance e la robustezza dell’algoritmo implementato.

## Struttura del Progetto

- **Codice sorgente**: Implementazione in C della tecnica KmerCO e dei Bloom Filter.
- **Dataset**: File FASTQ e CSV per test e validazione.
- **Risultati**: Output testuali e grafici per l’analisi delle prestazioni e dell’accuratezza.

## Compilazione ed esecuzione del file Balaenoptera.fastq

1. **Compilazione**  
   Utilizza il `makefile` implementato sulla base dei file esistenti e utilizzato per compilare il progetto. Di seguito, vengono riportati i dettagli della compilazione:
   ```
   make clean
   make
   ```
2. **Esecuzione**  
   Esegui il programma principale passando il file FASTQ desiderato (consigliato `balaenoptera.fastq` per test completi):
   ```
   ./KmerCo.exe balaenoptera.fastq
   ```
3. **Analisi dei risultati**  
   I risultati verranno salvati nei file di output (`Result.txt`, `Distinct.txt`, ecc.) e potranno essere visualizzati tramite gli script Python e i notebook presenti nella cartella `graphs/`.