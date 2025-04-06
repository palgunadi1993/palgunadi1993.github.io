document.addEventListener('DOMContentLoaded', function() {
  // Dark mode toggle
  const themeToggle = document.querySelector('.theme-toggle');
  const themeIcon = themeToggle.querySelector('i');
  
  // Check for saved theme preference or use preferred color scheme
  const savedTheme = localStorage.getItem('theme');
  const prefersDark = window.matchMedia('(prefers-color-scheme: dark)').matches;
  
  if (savedTheme === 'dark' || (!savedTheme && prefersDark)) {
    document.body.classList.add('dark-theme');
    themeIcon.classList.replace('fa-moon', 'fa-sun');
  }
  
  // Toggle theme
  themeToggle.addEventListener('click', function() {
    document.body.classList.toggle('dark-theme');
    
    if (document.body.classList.contains('dark-theme')) {
      localStorage.setItem('theme', 'dark');
      themeIcon.classList.replace('fa-moon', 'fa-sun');
    } else {
      localStorage.setItem('theme', 'light');
      themeIcon.classList.replace('fa-sun', 'fa-moon');
    }
  });
  
  // Highlight active page in navigation
  const currentPage = window.location.pathname.split("/").pop();
  const navLinks = document.querySelectorAll('nav a');
  
  navLinks.forEach(link => {
    const linkPage = link.getAttribute('href');
    if (linkPage === currentPage || (currentPage === '' && linkPage === 'index.html')) {
      link.classList.add('active');
    }
  });
  
  // Read More functionality for news items
  const readMoreLinks = document.querySelectorAll('.read-more-link');
  
  readMoreLinks.forEach(link => {
    link.addEventListener('click', function(e) {
      e.preventDefault();
      
      const newsItem = this.closest('.news-content');
      const details = newsItem.querySelector('.news-details');
      
      if (details.style.display === 'block') {
        details.style.display = 'none';
        this.innerHTML = 'Read more <i class="fas fa-chevron-right"></i>';
      } else {
        details.style.display = 'block';
        this.innerHTML = 'Read less <i class="fas fa-chevron-up"></i>';
      }
    });
  });

  // Toggle graphical abstracts
  const abstractToggles = document.querySelectorAll('.abstract-toggle');
  
  abstractToggles.forEach(toggle => {
    toggle.addEventListener('click', function() {
      const pubItem = this.closest('.publication-item');
      const abstract = pubItem.querySelector('.graphical-abstract');
      
      if (abstract.classList.contains('hidden-abstract')) {
        abstract.classList.remove('hidden-abstract');
        this.innerHTML = 'Hide graphical abstract <i class="fas fa-chevron-up"></i>';
      } else {
        abstract.classList.add('hidden-abstract');
        this.innerHTML = 'Show graphical abstract <i class="fas fa-chevron-down"></i>';
      }
    });
  });
  
  // Handle contact form submission if it exists
  const contactForm = document.querySelector('.contact-form');
  if (contactForm) {
    contactForm.addEventListener('submit', function(e) {
      e.preventDefault();
      // You can add form submission logic here
      // For now, we'll just show an alert
      alert('Form submission is not set up yet. Please configure your form handling.');
    });
  }
});
