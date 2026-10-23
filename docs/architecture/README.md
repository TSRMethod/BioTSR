# Architecture decisions

[`DECISIONS.md`](DECISIONS.md) holds BioTSR's architecture decision records,
ADR-001 to ADR-114. Each one records a single decision and its scientific or
engineering rationale. Together they document **why BioTSR 1.0 behaves as it
does**: why `math.pi` is the angle convention, why a key space is a
fingerprint, why SSE changes no key, why the command line holds no scientific
default, why each method runs under one standard profile, and so on.

For a present-tense overview of the package, read
[`../ARCHITECTURE.md`](../ARCHITECTURE.md) first. Come here when you want the
reason behind a specific choice.

## Reading the records

* The numbering is stable: source comments and tests cite these numbers, so
  every number resolves, including records that no longer describe current
  behaviour.
* **Accepted** records are in force. **Accepted, amended for 1.0** records are in
  force and have been restated for the 1.0 package.
* A record **Superseded by ADR-xxx** was replaced by the record named.
* A record marked **Historical: compatibility functionality removed in 1.0**
  describes a convention of the historical TSR programs. BioTSR 1.0 implements
  standard conventions only (ADR-087); the record says what the historical
  convention was and why the standard differs.
* A changed decision gets a new record rather than a silent edit.
