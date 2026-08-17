-- Step 00: create a dataset in your own project to hold intermediate tables.
-- Run once. Must be US multi-region to join against nih-sra-datastore.
--
-- Sandbox note: tables created here expire automatically after 60 days.
-- Steps 01 and 02 are cheap to re-run if that happens.

CREATE SCHEMA IF NOT EXISTS sra_ad OPTIONS(location = 'US');
