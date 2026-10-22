document.addEventListener("DOMContentLoaded", () => {

    const downloadButton =
        document.getElementById("download-config");

    downloadButton.addEventListener(
        "click",
        generateConfig
    );

});

(function () {
        const container = document.getElementById('LibraryIdentifiersContainer');
        const addBtn = document.getElementById('addLibraryIdentifier');
        const serialized = document.getElementById('LibraryIdentifiersSerialized');
        const warning = document.getElementById('libraryIdentifierWarning');
        const ngsfilesContainer = document.getElementById('ngsfilterFilesContainer');
        const multilineContainer = document.getElementById('multilineLibrariesContainer');
        const multilineToggle = document.getElementById('multilineLibraries');
        const r1Container = document.getElementById('r1FilesByLibraryContainer');
        const r1Serialized = document.getElementById('r1FilesSerialized');
        const r2Container = document.getElementById('r2FilesByLibraryContainer');
        const r2Serialized = document.getElementById('r2FilesSerialized');

        if (!container || !addBtn || !serialized) {
            return;
        }

        let r1State = {};
        let r2State = {};

        function getLibraryValues() {
            return Array.from(container.querySelectorAll('.LibraryIdentifier'))
                .map(input => input.value.trim())
                .filter(Boolean);
        }

        function getRequiredRowCount(libraries) {
            if (libraries.length === 1 && multilineToggle && multilineToggle.checked) {
                return 2;
            }

            return 1;
        }

        function normalizeStateForLibraries(state, libraries, rowCount, suffix) {
            const nextState = {};

            libraries.forEach(libraryId => {
                const previous = state[libraryId] || [];
                const values = [];

                for (let index = 0; index < rowCount; index += 1) {
                    values.push(previous[index] || `${libraryId}_${suffix}.fastq.gz`);
                }

                nextState[libraryId] = values;
            });

            return nextState;
        }

        function serializeFileState(state, libraries) {
            const result = [];

            libraries.forEach(libraryId => {
                const values = state[libraryId] || [];

                values.forEach(fileName => {
                    result.push({
                        library: libraryId,
                        file: fileName.trim()
                    });
                });
            });

            return result;
        }

        function renderFiles(containerEl, libraries, state, cssClass, suffix, rowCount) {
            if (!containerEl) {
                return;
            }

            containerEl.innerHTML = '';

            if (libraries.length === 0) {
                const helper = document.createElement('p');
                helper.textContent = '*Add a library identifier to generate file names.';
                containerEl.appendChild(helper);
                return;
            }

            libraries.forEach(libraryId => {
                for (let index = 0; index < rowCount; index += 1) {
                    const row = document.createElement('div');
                    row.style.display = 'grid';
                    row.style.gap = '4px';
                    row.style.marginBottom = '8px';

                    const label = document.createElement('label');
                    label.textContent = `${suffix} file for ${libraryId}`;

                    const input = document.createElement('input');
                    input.type = 'text';
                    input.className = cssClass;
                    input.dataset.library = libraryId;
                    input.dataset.index = String(index);
                    input.value = (state[libraryId] && state[libraryId][index]) || `${libraryId}_${suffix}.fastq.gz`;

                    row.appendChild(label);
                    row.appendChild(input);
                    containerEl.appendChild(row);
                }
            });
        }

        function bindFileInputs(containerEl, state, libraries, syncFn) {
            if (!containerEl) {
                return;
            }

            containerEl.querySelectorAll('input[type="text"]').forEach(input => {
                input.addEventListener('input', function () {
                    const libraryId = input.dataset.library;
                    const valueIndex = Number(input.dataset.index);

                    if (!state[libraryId]) {
                        state[libraryId] = [];
                    }

                    state[libraryId][valueIndex] = input.value;
                    syncFn(libraries);
                });
            });
        }

        function renderNgsfiles(libraries) {
            if (!ngsfilesContainer) {
                return;
            }

            ngsfilesContainer.innerHTML = '';

            if (libraries.length === 0) {
                const helper = document.createElement('p');
                helper.textContent = '*Add a library identifier to generate ngsfile names.';
                ngsfilesContainer.appendChild(helper);
                return;
            }

            libraries.forEach(libraryId => {
                const fileName = document.createElement('h5');
                fileName.textContent = `ngsfile_${libraryId}.tsv`;
                ngsfilesContainer.appendChild(fileName);
            });

            const helper = document.createElement('p');
            helper.textContent = '*Prepare one file for each library and rename them according to the library identifier.';
            ngsfilesContainer.appendChild(helper);
        }

        function validateLibraryIdentifiers() {
            const inputs = Array.from(container.querySelectorAll('.LibraryIdentifier'));
            const invalidValues = [];

            inputs.forEach(input => {
                const value = input.value.trim();
                const isValid = value.length === 0 || value.length === 4;

                if (!isValid) {
                    invalidValues.push(value);
                }

                input.setCustomValidity(
                    isValid ? '' : 'Library identifiers must have exactly 4 characters.'
                );
            });

            if (!warning) {
                return invalidValues.length === 0;
            }

            if (invalidValues.length > 0) {
                warning.style.display = 'block';
                warning.textContent = 'Library identifiers must have exactly 4 characters.';
                return false;
            }

            warning.style.display = 'none';
            warning.textContent = '';
            return true;
        }

        function syncLibraryIdentifiers(librariesParam) {
            const libraries = librariesParam || getLibraryValues();
            const rowCount = getRequiredRowCount(libraries);

            serialized.value = JSON.stringify(libraries);
            validateLibraryIdentifiers();
            renderNgsfiles(libraries);

            r1State = normalizeStateForLibraries(r1State, libraries, rowCount, 'R1');
            r2State = normalizeStateForLibraries(r2State, libraries, rowCount, 'R2');

            renderFiles(r1Container, libraries, r1State, 'r1-library-file', 'R1', rowCount);
            renderFiles(r2Container, libraries, r2State, 'r2-library-file', 'R2', rowCount);

            bindFileInputs(r1Container, r1State, libraries, syncLibraryIdentifiers);
            bindFileInputs(r2Container, r2State, libraries, syncLibraryIdentifiers);

            if (r1Serialized) {
                r1Serialized.value = JSON.stringify(serializeFileState(r1State, libraries));
            }

            if (r2Serialized) {
                r2Serialized.value = JSON.stringify(serializeFileState(r2State, libraries));
            }

            if (multilineContainer && multilineToggle) {
                if (libraries.length === 1) {
                    multilineContainer.style.display = 'block';
                } else {
                    multilineToggle.checked = false;
                    multilineContainer.style.display = 'none';
                }
            }
        }

        function buildRow(value = '') {
            const row = document.createElement('div');
            row.style.display = 'flex';
            row.style.gap = '8px';
            row.style.marginBottom = '6px';

            const input = document.createElement('input');
            input.type = 'text';
            input.className = 'LibraryIdentifier';
            input.placeholder = 'LIB2';
            input.value = value;

            const removeBtn = document.createElement('button');
            removeBtn.type = 'button';
            removeBtn.textContent = '−';
            removeBtn.setAttribute('aria-label', 'Remove this Library identifier');
            removeBtn.addEventListener('click', function () {
                row.remove();
                syncLibraryIdentifiers();
            });

            input.addEventListener('input', syncLibraryIdentifiers);

            row.appendChild(input);
            row.appendChild(removeBtn);
            return row;
        }

        addBtn.addEventListener('click', function () {
            container.appendChild(buildRow(''));
            syncLibraryIdentifiers();
        });

        container.querySelectorAll('.LibraryIdentifier').forEach(input => {
            input.addEventListener('input', syncLibraryIdentifiers);
        });

        if (multilineToggle) {
            multilineToggle.addEventListener('change', syncLibraryIdentifiers);
        }

        syncLibraryIdentifiers();
    })();

function generateConfig() {

    const libraryIdentifiersContainer =
        document.getElementById("LibraryIdentifiersContainer");

    if (libraryIdentifiersContainer) {
        const libraryInputs =
            Array.from(libraryIdentifiersContainer.querySelectorAll('.LibraryIdentifier'));

        const invalidInput =
            libraryInputs.find(input => {
                const value = input.value.trim();
                return value.length > 0 && value.length !== 4;
            });

        if (invalidInput) {
            invalidInput.reportValidity();
            invalidInput.focus();
            return;
        }
    }

    const project =
        document.getElementById("projectName").value;

    const r1Files =
        document.getElementById("r1FilesSerialized").value;

    const r2Files =
        document.getElementById("r2FilesSerialized").value;

    const threads =
        document.getElementById("threads").value;

    const quality =
        document.getElementById("quality").value;

    const config = `
# MJOLNIR3 Configuration

PROJECT=${project}

R1_FILES=${r1Files}

R2_FILES=${r2Files}

THREADS=${threads}

QUALITY=${quality}
`;

    const blob =
        new Blob([config], {
            type: "text/plain"
        });

    const link =
        document.createElement("a");

    link.href =
        URL.createObjectURL(blob);

    link.download =
        "mjolnir3_config.txt";

    link.click();

    URL.revokeObjectURL(link.href);
}