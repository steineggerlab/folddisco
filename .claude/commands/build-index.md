Guide the user through building a folddisco index and deploying it to the mmseqs2-app.

## Prerequisites

The folddisco databases live on an NFS disk. Mount it to your VM first:

```bash
sudo mkdir -p /mnt/folddisco-data
sudo mount -t nfs 10.165.165.250:/shared/folddisco-data /mnt/folddisco-data
```

Verify the mount:
```bash
ls /mnt/folddisco-data/
```

You should see existing databases and structure directories.

## 1. Prepare structure files

Place PDB or CIF files in a directory under `/mnt/folddisco-data/`. For example:

```
/mnt/folddisco-data/my_structures/
```

## 2. Build the folddisco binary

From the folddisco repo root:

```bash
cargo build --release
```

The `foldcomp` feature is enabled by default. The binary lands at `target/release/folddisco`.

## 3. Build the index

Run the index command pointing at the structure directory. Use **absolute paths** to avoid the
CWD-relative path bug (folddisco stores paths from `load_path()` verbatim in the `.lookup` file,
and worker pods run with a different CWD than where you built the index):

```bash
./target/release/folddisco index \
  -d /mnt/folddisco-data/my_structures/ \
  -o /mnt/folddisco-data/my_structures_folddisco \
  -t 16 \
  --verbose
```

This produces:
- `my_structures_folddisco.offset` — the hash index
- `my_structures_folddisco.lookup` — maps numeric IDs to file paths
- `my_structures_folddisco.type` — index config (hash type, bin counts, input format, etc.)

## 4. Verify lookup paths are absolute

**This is critical.** The worker pods mount the NFS at `/data/`, so lookup paths must be absolute
and rooted at `/data/`. Inspect the lookup file:

```bash
head -5 /mnt/folddisco-data/my_structures_folddisco.lookup
```

If paths are relative (e.g. `my_structures/foo.cif` instead of `/data/my_structures/foo.cif`),
fix them with the provided script:

```bash
python scripts/fix_lookup_paths.py /mnt/folddisco-data/my_structures_folddisco.lookup /data/
```

Dry-run first with `--dry-run` to preview changes.

## 5. Create the .params file

The mmseqs2-app API needs a `.params` JSON file to register the database:

```json
{
  "name": "MY_STRUCTURES",
  "version": "",
  "path": "my_structures_folddisco",
  "default": false,
  "order": 1,
  "taxonomy": false,
  "complex": false,
  "motif": false,
  "full_header": false,
  "index": "",
  "search": "",
  "multimer": "",
  "interface": false,
  "status": "COMPLETE"
}
```

Write it to the NFS mount:

```bash
cat > /mnt/folddisco-data/my_structures_folddisco.params << 'EOF'
{ ... }
EOF
```

## 6. Restart the API pod

```bash
kubectl rollout restart deployment/folddisco-api -n folddisco
```

Verify the new database appears:

```bash
curl http://folddisco-web.playground.genesistherapeutics.ai/api/databases
```

## Architecture notes

- **NFS mount on VM:** `10.165.165.250:/shared/folddisco-data` -> `/mnt/folddisco-data`
- **NFS mount in pods:** same NFS share mounted at `/data`
- **The path bug:** `load_path()` stores paths relative to CWD. Worker pods run with CWD `/`,
  so relative paths like `my_structures/foo.cif` resolve to `/my_structures/foo.cif` (wrong).
  Always use absolute paths or fix the lookup after building.
- **FCZDB indexes don't have this problem** because they read structures from the binary foldcomp
  DB file directly, not from individual file paths in the lookup.
