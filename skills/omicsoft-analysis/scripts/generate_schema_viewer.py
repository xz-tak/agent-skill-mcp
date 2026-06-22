#!/usr/bin/env python3
"""
Generate Standalone Schema Viewer
Creates a single HTML file with embedded JSON data that works without a web server

Usage:
    python generate_schema_viewer.py --json schema_report.json --output schema_viewer_standalone.html
"""

import argparse
import json
import sys


def generate_html(json_data):
    """Generate HTML with embedded JSON data"""

    # Escape and embed the JSON data
    json_str = json.dumps(json_data).replace('</script>', '<\\/script>')

    data_type = json_data.get('file_info', {}).get('data_type', 'deg')
    if data_type == 'expr':
        page_title = 'Expression Schema Browser'
        page_subtitle = 'Comprehensive exploration of Omicsoft expression data metadata'
    elif data_type == 'deg':
        page_title = 'DEG Schema Browser'
        page_subtitle = 'Comprehensive exploration of Omicsoft DEG analysis metadata'
    else:
        page_title = 'H5AD Schema Browser'
        page_subtitle = 'Comprehensive exploration of Omicsoft h5ad file metadata'

    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{page_title}</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: #f5f5f5;
            padding: 20px;
            line-height: 1.6;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            overflow: hidden;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
        }}
        .header h1 {{ font-size: 28px; margin-bottom: 10px; }}
        .header p {{ opacity: 0.9; font-size: 14px; }}
        .controls {{
            padding: 20px;
            background: #f8f9fa;
            border-bottom: 1px solid #dee2e6;
            display: flex;
            gap: 15px;
            flex-wrap: wrap;
            align-items: center;
        }}
        .search-box {{ flex: 1; min-width: 300px; }}
        .search-box input {{
            width: 100%;
            padding: 10px 15px;
            border: 1px solid #ced4da;
            border-radius: 4px;
            font-size: 14px;
        }}
        .filter-select {{
            padding: 10px 15px;
            border: 1px solid #ced4da;
            border-radius: 4px;
            font-size: 14px;
            background: white;
            cursor: pointer;
        }}
        .content {{ padding: 20px; }}
        .info-card {{
            background: #e3f2fd;
            border-left: 4px solid #2196f3;
            padding: 15px;
            margin-bottom: 20px;
            border-radius: 4px;
        }}
        .info-card h3 {{ color: #1976d2; margin-bottom: 10px; font-size: 16px; }}
        .info-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 10px;
        }}
        .info-item {{
            display: flex;
            justify-content: space-between;
            font-size: 14px;
        }}
        .category-section {{ margin-bottom: 30px; }}
        .category-header {{
            background: #495057;
            color: white;
            padding: 12px 15px;
            cursor: pointer;
            border-radius: 4px;
            display: flex;
            justify-content: space-between;
            align-items: center;
            transition: background 0.2s;
        }}
        .category-header:hover {{ background: #343a40; }}
        .category-header h2 {{ font-size: 18px; font-weight: 600; }}
        .category-badge {{
            background: rgba(255,255,255,0.2);
            padding: 4px 12px;
            border-radius: 12px;
            font-size: 12px;
        }}
        .category-content {{
            display: none;
            padding: 15px;
            background: #f8f9fa;
            border: 1px solid #dee2e6;
            border-top: none;
            border-radius: 0 0 4px 4px;
        }}
        .category-content.active {{ display: block; }}
        .column-item {{
            background: white;
            border: 1px solid #dee2e6;
            border-radius: 4px;
            margin-bottom: 15px;
            overflow: hidden;
        }}
        .column-header {{
            padding: 12px 15px;
            background: #f8f9fa;
            cursor: pointer;
            display: flex;
            justify-content: space-between;
            align-items: center;
            border-bottom: 1px solid #dee2e6;
        }}
        .column-header:hover {{ background: #e9ecef; }}
        .column-name {{ font-weight: 600; color: #495057; font-size: 15px; }}
        .column-stats {{
            display: flex;
            gap: 15px;
            font-size: 12px;
            color: #6c757d;
        }}
        .stat-badge {{
            background: #e9ecef;
            padding: 3px 8px;
            border-radius: 3px;
        }}
        .column-content {{
            display: none;
            padding: 15px;
            max-height: 500px;
            overflow-y: auto;
        }}
        .column-content.active {{ display: block; }}
        .values-table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 13px;
        }}
        .values-table thead {{
            background: #f8f9fa;
            position: sticky;
            top: 0;
        }}
        .values-table th {{
            padding: 10px;
            text-align: left;
            font-weight: 600;
            color: #495057;
            border-bottom: 2px solid #dee2e6;
        }}
        .values-table td {{
            padding: 8px 10px;
            border-bottom: 1px solid #f1f3f5;
        }}
        .values-table tbody tr:hover {{ background: #f8f9fa; }}
        .value-cell {{ font-family: 'Courier New', monospace; color: #212529; }}
        .count-cell {{ color: #495057; text-align: right; }}
        .percentage-cell {{ color: #6c757d; text-align: right; }}
        .progress-bar {{
            height: 4px;
            background: #e9ecef;
            border-radius: 2px;
            overflow: hidden;
            margin-top: 4px;
        }}
        .progress-fill {{
            height: 100%;
            background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
            transition: width 0.3s;
        }}
        .toggle-icon {{ transition: transform 0.2s; }}
        .toggle-icon.open {{ transform: rotate(180deg); }}
        .highlight {{ background-color: yellow; font-weight: bold; }}
        .no-results {{ padding: 40px; text-align: center; color: #6c757d; }}
        .demographic-tag {{
            display: inline-block;
            background: #28a745;
            color: white;
            padding: 2px 8px;
            border-radius: 3px;
            font-size: 11px;
            margin-left: 8px;
        }}
        .key-filter-tag {{
            display: inline-block;
            background: #007bff;
            color: white;
            padding: 2px 8px;
            border-radius: 3px;
            font-size: 11px;
            margin-left: 8px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>{page_title}</h1>
            <p>{page_subtitle}</p>
        </div>

        <div class="controls">
            <div class="search-box">
                <input type="text" id="searchInput" placeholder="Search for columns or values...">
            </div>
            <select id="categoryFilter" class="filter-select">
                <option value="all">All Categories</option>
                <option value="key_filtering">Key Filtering</option>
                <option value="demographic">Demographics</option>
                <option value="case_disease">Case Disease</option>
                <option value="control_disease">Control Disease</option>
                <option value="treatment">Treatment</option>
                <option value="response">Response</option>
                <option value="sample">Sample</option>
                <option value="comparison_details">Comparison Details</option>
                <option value="other">Other</option>
            </select>
            <button id="expandAll" class="filter-select">Expand All</button>
            <button id="collapseAll" class="filter-select">Collapse All</button>
        </div>

        <div class="content">
            <div id="fileInfo"></div>
            <div id="schemaContent"></div>
        </div>
    </div>

    <script>
        // Embedded schema data
        const schemaData = {json_str};
        let currentFilter = 'all';
        let searchTerm = '';

        // Initialize on load
        renderFileInfo();
        renderSchema();

        function renderFileInfo() {{
            const info = schemaData.file_info;
            const nVars = info.n_variables || (schemaData.var_info && schemaData.var_info.n_genes) || 'N/A';
            const nVarsStr = typeof nVars === 'number' ? nVars.toLocaleString() : nVars;
            const dataType = info.data_type || '';
            const typeBadge = dataType ? `<span style="background:#667eea;color:white;padding:2px 8px;border-radius:3px;font-size:11px;margin-left:8px;">${{dataType.toUpperCase()}}</span>` : '';
            const html = `
                <div class="info-card">
                    <h3>File Information ${{typeBadge}}</h3>
                    <div class="info-grid">
                        <div class="info-item">
                            <span>Observations:</span>
                            <strong>${{info.n_observations.toLocaleString()}}</strong>
                        </div>
                        <div class="info-item">
                            <span>Variables (Genes):</span>
                            <strong>${{nVarsStr}}</strong>
                        </div>
                    </div>
                    <p style="margin-top: 10px; font-size: 13px; color: #495057;">
                        ${{info.observation_description || ''}}
                    </p>
                </div>
            `;
            document.getElementById('fileInfo').innerHTML = html;
        }}

        function renderSchema() {{
            const categories = schemaData.column_categories;
            const knownNames = {{
                'key_filtering': 'Key Filtering Columns',
                'demographic': 'Demographics',
                'case_disease': 'Case Disease Information',
                'control_disease': 'Control Disease Information',
                'treatment': 'Treatment Information',
                'response': 'Response Information',
                'sample': 'Sample Information',
                'comparison_details': 'Comparison Details',
                'experiment': 'Experiment Information',
                'study': 'Study Information',
                'biological': 'Biological Annotations',
                'disease': 'Disease Information'
            }};
            const categoryNames = {{}};
            for (const key of Object.keys(categories)) {{
                categoryNames[key] = knownNames[key] || key.replace(/_/g, ' ').replace(/\\b\\w/g, c => c.toUpperCase());
            }}

            let html = '';

            for (const [catKey, columns] of Object.entries(categories)) {{
                if (columns.length === 0) continue;

                const categoryName = categoryNames[catKey] || catKey;
                const filteredColumns = filterColumns(columns, catKey);

                if (filteredColumns.length === 0 && currentFilter !== 'all') continue;

                html += `
                    <div class="category-section" data-category="${{catKey}}">
                        <div class="category-header" onclick="toggleCategory('${{catKey}}')">
                            <h2>${{categoryName}}</h2>
                            <div>
                                <span class="category-badge">${{filteredColumns.length}} columns</span>
                                <span class="toggle-icon" id="icon-${{catKey}}">▼</span>
                            </div>
                        </div>
                        <div class="category-content" id="content-${{catKey}}">
                            ${{renderColumns(filteredColumns, catKey)}}
                        </div>
                    </div>
                `;
            }}

            const allCategorizedCols = new Set(Object.values(categories).flat());
            const otherColumns = Object.keys(schemaData.obs_columns)
                .filter(col => !allCategorizedCols.has(col));

            if (otherColumns.length > 0 && (currentFilter === 'all' || currentFilter === 'other')) {{
                const filteredOther = filterColumns(otherColumns, 'other');
                if (filteredOther.length > 0) {{
                    html += `
                        <div class="category-section" data-category="other">
                            <div class="category-header" onclick="toggleCategory('other')">
                                <h2>Other Columns</h2>
                                <div>
                                    <span class="category-badge">${{filteredOther.length}} columns</span>
                                    <span class="toggle-icon" id="icon-other">▼</span>
                                </div>
                            </div>
                            <div class="category-content" id="content-other">
                                ${{renderColumns(filteredOther, 'other')}}
                            </div>
                        </div>
                    `;
                }}
            }}

            document.getElementById('schemaContent').innerHTML = html ||
                '<div class="no-results">No columns match your search criteria.</div>';
        }}

        function filterColumns(columns, category) {{
            if (currentFilter !== 'all' && currentFilter !== category) {{
                return [];
            }}

            if (!searchTerm) {{
                return columns;
            }}

            return columns.filter(col => {{
                const colData = schemaData.obs_columns[col];
                const colMatch = col.toLowerCase().includes(searchTerm.toLowerCase());
                const valueMatch = colData.all_values.some(v =>
                    v.value.toLowerCase().includes(searchTerm.toLowerCase())
                );
                return colMatch || valueMatch;
            }});
        }}

        function renderColumns(columns, category) {{
            return columns.map(col => {{
                const colData = schemaData.obs_columns[col];
                const isDemographic = category === 'demographic';
                const isKeyFilter = category === 'key_filtering';

                return `
                    <div class="column-item">
                        <div class="column-header" onclick="toggleColumn('${{col}}')">
                            <div>
                                <span class="column-name">${{col}}</span>
                                ${{isDemographic ? '<span class="demographic-tag">DEMO</span>' : ''}}
                                ${{isKeyFilter ? '<span class="key-filter-tag">FILTER</span>' : ''}}
                            </div>
                            <div class="column-stats">
                                <span class="stat-badge">Unique: ${{colData.unique_values}}</span>
                                <span class="stat-badge">Total: ${{colData.total_values}}</span>
                                <span class="toggle-icon" id="icon-col-${{col}}">▼</span>
                            </div>
                        </div>
                        <div class="column-content" id="content-col-${{col}}">
                            ${{renderValues(colData.all_values)}}
                        </div>
                    </div>
                `;
            }}).join('');
        }}

        function renderValues(values) {{
            if (values.length === 0) {{
                return '<p style="padding: 10px; color: #6c757d;">No values</p>';
            }}

            let filteredValues = values;
            if (searchTerm) {{
                filteredValues = values.filter(v =>
                    v.value.toLowerCase().includes(searchTerm.toLowerCase())
                );
            }}

            if (filteredValues.length === 0) {{
                return '<p style="padding: 10px; color: #6c757d;">No matching values</p>';
            }}

            const rows = filteredValues.map(v => {{
                const valueText = highlightSearch(v.value);
                return `
                    <tr>
                        <td class="value-cell">${{valueText}}</td>
                        <td class="count-cell">${{v.count.toLocaleString()}}</td>
                        <td class="percentage-cell">
                            ${{v.percentage.toFixed(2)}}%
                            <div class="progress-bar">
                                <div class="progress-fill" style="width: ${{Math.min(v.percentage, 100)}}%"></div>
                            </div>
                        </td>
                    </tr>
                `;
            }}).join('');

            return `
                <table class="values-table">
                    <thead>
                        <tr>
                            <th>Value</th>
                            <th style="text-align: right;">Count</th>
                            <th style="text-align: right;">Percentage</th>
                        </tr>
                    </thead>
                    <tbody>
                        ${{rows}}
                    </tbody>
                </table>
            `;
        }}

        function highlightSearch(text) {{
            if (!searchTerm) return text;
            const regex = new RegExp(`(${{searchTerm}})`, 'gi');
            return text.replace(regex, '<span class="highlight">$1</span>');
        }}

        function toggleCategory(category) {{
            const content = document.getElementById(`content-${{category}}`);
            const icon = document.getElementById(`icon-${{category}}`);

            if (content.classList.contains('active')) {{
                content.classList.remove('active');
                icon.classList.remove('open');
            }} else {{
                content.classList.add('active');
                icon.classList.add('open');
            }}
        }}

        function toggleColumn(column) {{
            const content = document.getElementById(`content-col-${{column}}`);
            const icon = document.getElementById(`icon-col-${{column}}`);

            if (content.classList.contains('active')) {{
                content.classList.remove('active');
                icon.classList.remove('open');
            }} else {{
                content.classList.add('active');
                icon.classList.add('open');
            }}
        }}

        document.getElementById('searchInput').addEventListener('input', (e) => {{
            searchTerm = e.target.value;
            renderSchema();
        }});

        document.getElementById('categoryFilter').addEventListener('change', (e) => {{
            currentFilter = e.target.value;
            renderSchema();
        }});

        document.getElementById('expandAll').addEventListener('click', () => {{
            document.querySelectorAll('.category-content').forEach(el => {{
                el.classList.add('active');
            }});
            document.querySelectorAll('.toggle-icon').forEach(el => {{
                el.classList.add('open');
            }});
        }});

        document.getElementById('collapseAll').addEventListener('click', () => {{
            document.querySelectorAll('.category-content, .column-content').forEach(el => {{
                el.classList.remove('active');
            }});
            document.querySelectorAll('.toggle-icon').forEach(el => {{
                el.classList.remove('open');
            }});
        }});
    </script>
</body>
</html>'''

    return html


def main():
    parser = argparse.ArgumentParser(
        description='Generate standalone schema viewer HTML with embedded JSON'
    )

    parser.add_argument(
        '--json',
        required=True,
        help='Path to schema_report.json'
    )

    parser.add_argument(
        '--output',
        required=True,
        help='Output HTML file path'
    )

    args = parser.parse_args()

    try:
        # Load JSON
        print(f"Loading {args.json}...")
        with open(args.json, 'r') as f:
            data = json.load(f)

        print(f"✓ Loaded {len(data['obs_columns'])} columns")

        # Generate HTML
        print(f"Generating HTML...")
        html = generate_html(data)

        # Write output
        with open(args.output, 'w') as f:
            f.write(html)

        print(f"✓ Generated {args.output} ({len(html):,} bytes)")
        print(f"\nYou can now open {args.output} directly in your browser!")

    except Exception as e:
        print(f"✗ Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
