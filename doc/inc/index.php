<?php include('faqdata.php'); ?>
<?php
//Output a table of contents with popping up

$sectionCounter = 1;
foreach ($faq as $section) {
  echo '<h4>'.$section['title'].'</h4>';
  echo '<ol>';
    foreach($section['content'] as $item) {
    $label = $item['label'];
  ?>
  <li>
    <div id="<?php echo $label ?>" class="reveal">
      <?php echo $item['question']; ?>
    </div>
    <div id="<?php echo $label ?>_ans" class="answer">
       <?php echo $item['answer'];?>
    </div>
  </li>
  <?php
  }
  echo "</ol><br/>\n";
}
?>
