#!/usr/bin/env python3
"""Export Mermaid and an nf-metro SVG from Nextflow's existing HTML DAG."""

import argparse
from html.parser import HTMLParser
from pathlib import Path
import shutil
import re
import subprocess
import sys


class MermaidDiagram(HTMLParser):
    def __init__(self):
        super().__init__()
        self.active = False
        self.parts = []

    def handle_starttag(self, tag, attrs):
        if tag == "pre" and "mermaid" in dict(attrs).get("class", "").split():
            self.active = True

    def handle_endtag(self, tag):
        if tag == "pre":
            self.active = False

    def handle_data(self, data):
        if self.active:
            self.parts.append(data)


def colour_sections(output):
    """Add section colours after rendering so routing and spacing stay intact."""
    palette = {
        "database_generation": ("#16877b", "#67cfba", "#e8f5f1", "#293e38"),
        "assessment": ("#b66b16", "#f1b65b", "#fff3df", "#43382b"),
        "statistics": ("#7250b0", "#b6a0ed", "#f1ebfa", "#383146"),
    }
    rules = []
    for section, (light, dark, background, dark_background) in palette.items():
        selector = f'[data-section-id="{section}"]'
        colour = f"light-dark({light}, {dark})"
        rules.extend([
            f'.nf-metro-section-box{selector} {{ fill: light-dark({background}, {dark_background}); stroke: {colour}; }}',
            f'.nf-metro-section-label{selector}, .nf-metro-section-num-circle{selector} {{ fill: {colour}; }}',
            f'.nf-metro-station{selector} {{ fill: {colour}; }}',
        ])
    svg = output.read_text()
    svg = svg.replace('</svg>', '<style id="virmarkdb-colours">\n' + '\n'.join(rules) + '\n</style>\n</svg>')
    output.write_text(svg)


def simplify_routes(text):
    """Hide transitive links in the overview; keep every station and reachability."""
    pattern = re.compile(r"^\s*(\w+) -->\|\w+\| (\w+)\s*$")
    edges = [match.groups() for line in text.splitlines()
             if (match := pattern.match(line))]
    adjacency = {}
    for source, target in edges:
        adjacency.setdefault(source, set()).add(target)

    def redundant(source, target):
        pending = list(adjacency[source] - {target})
        seen = set()
        while pending:
            node = pending.pop()
            if node == target:
                return True
            if node not in seen:
                seen.add(node)
                pending.extend(adjacency.get(node, ()))
        return False

    return "\n".join(line for line in text.splitlines()
                     if not ((match := pattern.match(line)) and redundant(*match.groups())))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dag", type=Path, help="Nextflow HTML or Mermaid DAG")
    parser.add_argument("--vertical", action="store_true", help="Stack sections vertically")
    args = parser.parse_args()
    dag = args.dag
    if dag.suffix == ".html":
        diagram = MermaidDiagram()
        diagram.feed(dag.read_text())
        if not diagram.parts:
            raise ValueError(f"No Mermaid diagram found in {dag}")
        source = dag.with_suffix(".mmd")
        source.write_text("".join(diagram.parts).strip() + "\n")
    elif dag.suffix == ".mmd":
        source = dag
    else:
        raise ValueError("Automatic metro rendering requires an HTML or MMD DAG")

    print(f"Mermaid diagram: {source}", flush=True)
    executable = shutil.which("nf-metro")
    if not executable:
        print("nf-metro is unavailable; HTML/Mermaid exports are still available. "
              "Install it with: pip install nf-metro", file=sys.stderr)
        return

    # nf-metro 3 uses --from nextflow; older releases use --from-nextflow.
    help_text = subprocess.run([executable, "render", "--help"], check=True,
                               capture_output=True, text=True).stdout
    source_option = ["--from-nextflow"] if "--from-nextflow" in help_text else ["--from", "nextflow"]
    output = source.with_name("virmarkdb-metro.svg")
    if args.vertical:
        styled = source.with_name("virmarkdb-metro-vertical.mmd")
        subprocess.run([executable, "convert", str(source), "--from", "nextflow",
                        "--title", "VirMarkDB", "-o", str(styled)], check=True)
        text = styled.read_text()
        sections = re.findall(r"^\s*subgraph (\w+) \[", text, re.MULTILINE)
        # Colour each dependency by its destination branch, preserving all edges.
        nodes = {}
        section = None
        lines = []
        for line in text.splitlines():
            if line.startswith("%%metro line:"):
                continue
            match = re.match(r"\s*subgraph (\w+) ", line)
            if match:
                section = match[1]
            station = re.match(r"\s*(\w+)\(\[", line)
            if station:
                nodes[station[1]] = section
            lines.append(line)
        text = re.sub(r"(\w+) -->\|\w+\| (\w+)",
                      lambda edge: f"{edge[1]} -->|{nodes[edge[2]]}| {edge[2]}",
                      "\n".join(lines))
        text = simplify_routes(text)
        notes = []
        # In this presentation view, show the early metadata dependency as a note
        # rather than a long route crossing the whole generation section.
        marker_link = "generate_manifest -->|assessment| collect_pool"
        if marker_link in text:
            text = "\n".join(line for line in text.splitlines() if line.strip() != marker_link)
            notes.append("%%metro caption: Collect Pool also receives marker definitions from Generate Manifest.")
        palette = {"database_generation": "#2db572", "assessment": "#eba642",
                   "statistics": "#a48ae0"}
        routes = [f"%%metro line: {section} | {section.replace('_', ' ').title()} | {palette.get(section, '#2db572')}"
                  for section in sections]
        if sections == ["database_generation", "assessment", "statistics"]:
            grid = ["%%metro grid: database_generation | 0,0,1,2",
                    "%%metro grid: assessment | 0,1", "%%metro grid: statistics | 1,1"]
        else:
            grid = [f"%%metro grid: {section} | 0,{row}" for row, section in enumerate(sections)]
        styled.write_text("\n".join(grid + routes + notes) + "\n" + text + "\n")
        subprocess.run([executable, "render", str(styled), "-o", str(output), "--validate"], check=True)
        colour_sections(output)
        return
    subprocess.run([executable, "render", str(source), "-o", str(output), *source_option], check=True)
    colour_sections(output)


if __name__ == "__main__":
    main()
