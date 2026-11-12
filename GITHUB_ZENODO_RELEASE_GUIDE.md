# Beginner upload instructions for the ACTINO–TransATor supporting dataset

1. Create an **empty** GitHub repository (suggested name: `ACTINO-TransATor-supplementary`) on your account. Choose *Private* initially if you need author review, and do not generate an extra README or template licence at creation.
2. Provide the repository link in your ChatGPT conversation. With the linked GitHub connector granted access, the prepared project files can be committed in their structured directories. Alternatively unzip the distribution locally and use GitHub Desktop/Git to commit all files; **uploading the ZIP as one ordinary GitHub file does not unpack it**.
3. Confirm file counts, six figures, twelve numbered supplementary tables, 18 comparative genomic neighborhoods and the executable script `reproduce_figures.py`. Do not upload your private MASTER v050 or any unpublished sensitive laboratory notes by default.
4. Review scientific records, authorship/credit, permissions and select a compatible licence before public release. GitHub provides licence templates. Do not create a CITATION.cff with guessed author names.
5. Change visibility to **Public** once author-approved; public visibility is needed for Zenodo/GitHub archiving.
6. Create or sign into Zenodo, link your GitHub account, go to the GitHub section, sync repositories and enable this repository **before** creating the release.
7. In GitHub choose Releases → Draft a new release → tag `v1.0.0` → publish. Check the Zenodo record and capture its genuine DOI.
8. Replace placeholders in `DATA_AVAILABILITY_TEXT_TEMPLATE.md` and the Word manuscript's Availability of data and materials section with the actual GitHub release URL and Zenodo DOI.
9. If the selected journal asks for Supplementary Files attached to the manuscript, submit its supplementary PDF and curated data ZIP separately as required. A GitHub link alone does not necessarily replace submission supplementary files.

Official guides:
- GitHub creating a repository: https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-new-repository
- GitHub adding files: https://docs.github.com/en/repositories/working-with-files/managing-files/adding-a-file-to-a-repository
- Zenodo enable repository: https://help.zenodo.org/docs/github/enable-repository/
- Zenodo archive GitHub release: https://help.zenodo.org/docs/github/archive-software/github-upload/
