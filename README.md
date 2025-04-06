# Kadek Hendrawan Palgunadi - Personal Academic Website

This repository contains the code for my personal academic website hosted at [palgunadi1993.github.io](https://palgunadi1993.github.io).

## Website Overview

This website showcases my academic profile as an Assistant Professor in Geophysical Engineering at Institut Teknologi Sepuluh Nopember and a Postdoctoral Researcher at ETH Zurich. It includes:

- Professional background and education
- Research interests and publications
- Details about the Earthquake Mechanics Research Group
- News and updates
- Useful links and resources
- Contact information

## Deployment Instructions

### Initial Setup

To deploy this website to GitHub Pages under the domain palgunadi1993.github.io:

1. **Create a GitHub repository** named exactly `palgunadi1993.github.io`
   - The repository name must match your GitHub username followed by .github.io

2. **Clone the repository** to your local machine
   ```bash
   git clone https://github.com/palgunadi1993/palgunadi1993.github.io.git
   cd palgunadi1993.github.io
   ```

3. **Copy all website files** to this directory
   - All HTML files should be in the root directory
   - Create subdirectories for assets: images, files, etc.

4. **Add, commit, and push** your files
   ```bash
   git add .
   git commit -m "Initial website deployment"
   git push origin main
   ```

5. **Enable GitHub Pages**
   - Go to your repository settings on GitHub
   - Navigate to the "Pages" section
   - Select the "main" branch as your source
   - Click Save

6. Wait a few minutes for GitHub to build and deploy your site. Your website will be available at `https://palgunadi1993.github.io`

### Adding Your Profile Photo

1. Create an `images` folder in the root directory
2. Add your photo to this folder (recommended format: JPG or PNG, size: 500x500px)
3. Update references in the HTML files:
   - `index.html`: Update `<img src="profile-placeholder.jpg"...>` to `<img src="images/your-photo.jpg"...>`
   - `group.html`: Also update the image path here

### Adding PDF Files

1. Create a `files` or `papers` folder in the root directory
2. Add your PDFs to this folder
3. Reference them in your HTML with relative paths: `<a href="files/your-paper.pdf"...>`

## Updating Content

### Adding News Items

1. Open `news.html` in a text editor
2. Add a new news item at the top of the news list:
   ```html
   <!-- Latest News Item -->
   <div class="news-item">
     <div class="news-date">
       <div class="date-box">
         <span class="month">Apr</span>
         <span class="year">2025</span>
       </div>
     </div>
     <div class="news-content">
       <h3>Your News Title</h3>
       <p>Brief description of your news...</p>
       <p class="read-more"><a href="#" class="read-more-link">Read more <i class="fas fa-chevron-right"></i></a></p>
       <div class="news-details">
         <p>Detailed content about your news. This will be hidden initially and shown when "Read more" is clicked.</p>
       </div>
     </div>
   </div>
   ```

### Adding Publications

1. Open `research.html` in a text editor
2. Add a new publication at the top of the publications list:
   ```html
   <!-- Latest Publication -->
   <div class="publication-item">
     <div class="pub-date">
       <div class="date-box">
         <span class="year">2025</span>
       </div>
     </div>
     <div class="pub-content">
       <h3 class="publication-title">Your Publication Title</h3>
       <p class="publication-authors">Author 1, Author 2, <strong>Kadek Hendrawan Palgunadi</strong>, Author 3</p>
       <p class="publication-venue">Journal Name, Vol(Issue), pp. XX-XX</p>
       <div class="publication-links">
         <a href="files/paper.pdf" target="_blank"><i class="fas fa-file-pdf"></i> PDF</a>
         <a href="https://github.com/your-repo" target="_blank"><i class="fas fa-code"></i> Code</a>
         <a href="https://doi.org/your-doi" target="_blank"><i class="fas fa-book"></i> DOI</a>
       </div>
     </div>
   </div>
   ```

### Modifying Website Content

1. Edit the HTML files directly for content changes
2. Test changes locally before pushing to GitHub:
   ```bash
   python -m http.server
   ```
   Then visit `http://localhost:8000` in your browser

### Updating Styles

1. Edit `styles.css` to change colors, fonts, spacing, etc.
2. The theme is controlled by CSS variables in the `:root` selector

## Regular Maintenance

To keep your website up to date:

1. Pull the latest changes (if working from multiple devices)
   ```bash
   git pull origin main
   ```

2. Make your changes to the appropriate HTML files

3. Test locally with a web server:
   ```bash
   python -m http.server
   ```

4. Commit and push your changes:
   ```bash
   git add .
   git commit -m "Description of your updates"
   git push origin main
   ```

5. GitHub Pages will automatically rebuild and deploy your updated site

## Troubleshooting

- If your site isn't updating after a push, check the "Actions" tab in your GitHub repository to see if there are any build errors
- Make sure all links and file paths are correct, especially for images and PDFs
- If using custom domains later on, you'll need to configure DNS settings and add a CNAME file

## Technologies Used

- HTML5
- CSS3
- JavaScript
- Font Awesome for icons
- GitHub Pages for hosting
